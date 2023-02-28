# %%
import rioxarray as rxr
import xarray as xr
from rasterio.plot import plotting_extent
from rasterio.transform import xy
import numpy as np
import geopandas as gpd
import cartopy.crs as ccrs
import cartopy
from pyproj import Transformer
import matplotlib.pyplot as plt
from collections import OrderedDict
from pathlib import Path
from shapely.geometry import Point, LineString
from itertools import pairwise

import warnings
warnings.simplefilter('ignore')

def generate_network(flow_dir_fn, stations_fn, DIST_PROJ, elevation_fn=None):
    fdr = rxr.open_rasterio(flow_dir_fn, masked=True)
    band = fdr.sel(band=1)

    band_vicfmt = band

    # create id of each cell
    band_ids = np.full_like(band_vicfmt, np.nan)
    _id = 0
    for i in np.arange(band_vicfmt.shape[0]):
        for j in np.arange(band_vicfmt.shape[1]):
            if not np.isnan(band_vicfmt[i][j]):
                band_ids[i][j] = _id
                _id += 1
    # create nodes
    # extract x-y points
    xs, _ = xy(fdr.rio.transform(), 0, np.arange(fdr.rio.width))
    _, ys = xy(fdr.rio.transform(), np.arange(fdr.rio.height), 9)

    XX, YY = np.meshgrid(xs, ys)

    # create gdf
    grid_pts = gpd.GeoDataFrame(data={'id': band_ids.flatten()}, geometry=[Point(x, y) for x, y in zip(XX.flatten(), YY.flatten())], crs='EPSG:4326')
    grid_pts['to_id'] = np.nan
    # create id of each cell
    operations = {
        1: [-1, 0],  # N
        2: [-1, 1],  # NE
        3: [0, 1],   # E
        4: [1, 1],   # SE
        5: [1, 0],   # S
        6: [1, -1],  # SW
        7: [0, -1],  # W
        8: [-1, -1], # NW
    }

    to_band_ids = np.full_like(band_vicfmt, np.nan)
    for x in np.arange(band_vicfmt.shape[0]):   
        for y in np.arange(band_vicfmt.shape[1]):  
            if not np.isnan(band_vicfmt[x][y]):
                if isinstance(band_vicfmt, xr.DataArray):
                    direction = int(band_vicfmt[x][y].values)
                else:
                    direction = band_vicfmt[x][y]
                op = operations[direction]
                new_x, new_y = x + op[0], y + op[1]
                if new_x < band_ids.shape[0] and new_y < band_ids.shape[1]:
                    # print(band_ids[x, y], (x, y), directions[direction], " > | >", band_ids[new_x, new_y], (new_x, new_y))
                    to_band_ids[x, y] = band_ids[new_x, new_y]
                else:
                    # print(band_ids[x, y], (x, y), directions[direction], " > | >", np.nan, (new_x, new_y))
                    to_band_ids[x, y] = np.nan

    grid_pts['to_id'] = to_band_ids.flatten()

    grid_pts.dropna(inplace=True)

    # Lets import the reservoir locations
    reservoirs = gpd.read_file(stations_fn)
    reservoirs['geometry'] = gpd.points_from_xy(reservoirs['lon'], reservoirs['lat'])
    reservoirs.set_crs('epsg:4326', inplace=True)

    reservoirs_to_join = reservoirs[['name', 'lon', 'lat', 'geometry']]

    if elevation_fn is not None:
        elevation = rxr.open_rasterio(elevation_fn)
        transformer = Transformer.from_crs("EPSG:4326", elevation.rio.crs, always_xy=True)
        indices = [transformer.transform(lon, lat) for (lon, lat) in zip(reservoirs['geometry'].x, reservoirs['geometry'].y)]
        reservoirs_to_join['elevation'] = [elevation.sel(x=xx, y=yy, method="nearest").values[0] for xx, yy in indices]

    grid_pts_joined_stns = grid_pts.to_crs(DIST_PROJ).sjoin_nearest(reservoirs_to_join.to_crs(DIST_PROJ), how='left', distance_col='dist [m]').to_crs('epsg:4326')
    # TODO: give a warning if multiple stations are within the same cell

    grid_pts_joined_stns['min_dist'] = grid_pts_joined_stns.groupby('name')['dist [m]'].transform('min')
    grid_pts_joined_stns.loc[grid_pts_joined_stns['dist [m]'] != grid_pts_joined_stns['min_dist'], ['index_right', 'name', 'lon', 'lat', 'dist [m]', 'min_dist', 'elevation']] = np.nan

    # save the root index
    root = None
    # Arrows
    for i, row in grid_pts.iterrows():
        from_loc = row['geometry']
        try:
            to_loc = grid_pts.loc[grid_pts['id'] == row['to_id']].iloc[0]['geometry']
        except IndexError:
            to_loc = from_loc
            root = row

        x, y = from_loc.x, from_loc.y
        dx, dy = to_loc.x - x, to_loc.y - y

    grid_pts_joined_stns = grid_pts_joined_stns.to_crs(DIST_PROJ)
    stations = grid_pts_joined_stns.dropna()

    dist = lambda left_pt, right_pt: np.sqrt((left_pt.x-right_pt.x)**2 + (left_pt.y-right_pt.y)**2)

    network = OrderedDict()
    for i, stn in stations.iterrows():
        id = stn['index_right']
        length = 0
        
        if id not in network:
            network[id] = OrderedDict([
                ('network', []),
                ('lengths', [])
            ])

        current_id = stn['id']
        while True:
            current_node = grid_pts_joined_stns.loc[grid_pts_joined_stns['id']==current_id]
            next_id = current_node['to_id'].iloc[0]
            next_node = grid_pts_joined_stns.loc[grid_pts_joined_stns['id']==next_id]
            # print("\n\n\n", current_node['geometry'], next_node['geometry'])
            if len(next_node) > 0:
                length += dist(current_node['geometry'].iloc[0], next_node['geometry'].iloc[0])
            if not np.isnan(current_node['index_right'].iloc[0]):
                network[id]['network'].append(current_node['index_right'].iloc[0])
                network[id]['lengths'].append(length)

            if current_id == root['id']:
                break

            current_id = next_id

    grid_pts_joined_stns = grid_pts_joined_stns.to_crs(DIST_PROJ)

    from_id = []
    to_id = []
    from_name = []
    to_name = []
    geoms = []
    lengths = []
    delta_elevs = []
    slopes = []

    for id in network:
        downstreams = network[id]['network']
        link_lengths = network[id]['lengths']

        for id_pair, length in zip(pairwise(downstreams), link_lengths[1:]):  # the first element of `link_lengths` is garbage value. It is the distance of a single link
            start_res = grid_pts_joined_stns.loc[grid_pts_joined_stns['index_right']==id_pair[0]].to_crs('epsg:4326').iloc[0]  # Get the points in EPSG:4326 CRS (lat-long)
            end_res = grid_pts_joined_stns.loc[grid_pts_joined_stns['index_right']==id_pair[1]].to_crs('epsg:4326').iloc[0]  # Get the points in EPSG:4326 CRS (lat-long)

            link = LineString([Point(pt.x, pt.y) for pt in (start_res.geometry, end_res.geometry)])

            geoms.append(link)
            from_id.append(id_pair[0])
            to_id.append(id_pair[1])
            from_name.append(start_res['name'])
            to_name.append(end_res['name'])
            lengths.append(length)
            delta_elevs.append(end_res['elevation']-start_res['elevation'])
            slopes.append((end_res['elevation']-start_res['elevation'])/length)

    res_network = gpd.GeoDataFrame(data={
        'from_id': from_id,
        'to_id': to_id,
        'from_name': from_name,
        'to_name': to_name,
        'geometry': geoms,
        'length': lengths,
        'delta_elevation': delta_elevs,
        'slope': slopes
    })

    # Save `res_network` for further use
    return res_network, grid_pts_joined_stns.dropna().to_crs('epsg:4326')


def main():
    DIST_PROJ = "+proj=eqdc +lon_0=-103.7988281 +lat_1=35.7127609 +lat_2=43.8942567 +lat_0=39.8035088 +datum=WGS84 +units=m +no_defs"
    flow_dir_fn = Path("data/colorado/basins/gunnison/ro/pars/fl.tif")
    stations_fn = Path("data/gunnison_reservoirs/gunnison_reservoirs_fdr_corrected_usgs_gages.csv")
    elevation_fn = Path("global_data/global_elevation_data/World_e-Atlas-UCSD_SRTM30-plus_v8.tif")
    save_fn = Path("data/gunnison_rivreg/rivreg_network.shp")
    save_fn.parent.mkdir(parents=True, exist_ok=True) # make parent directory if it doesn't exist

    assert flow_dir_fn.exists()
    assert stations_fn.exists()
    assert elevation_fn.exists()

    network, network_pts = generate_network(flow_dir_fn, stations_fn, DIST_PROJ, elevation_fn)

    network.to_file(save_fn)
    network_pts.to_file(save_fn.with_name('rivreg_network_pts.shp'))


if __name__ == '__main__':
    main()