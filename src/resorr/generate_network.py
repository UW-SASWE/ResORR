from pathlib import Path
import rioxarray as rxr
import xarray as xr
import numpy as np
import geopandas as gpd
import networkx as nx
import geonetworkx as gnx


import warnings
warnings.simplefilter('ignore')

def generate_network(
        flow_dir_fn, 
        stations_fn, 
        save_dir, 
        dist_proj=None, 
        elevation_fn=None
    ) -> (gpd.GeoDataFrame, gpd.GeoDataFrame): 
    """generate reservoir network using flow direction file and reservoir locations.

    Args:
        flow_dir_fn (str or Path): path to flow direction file
        stations_fn (str or Path): path to reservoir locations file
        save_dir (str or Path): path to save directory
        dist_proj (str, optional): projection to use for optionally calculating distance between reservoirs. Defaults to None.
        elevtion_fn (str or Path, optional): path to elevation file used for optionally adding elevation data to each reservoir. Defaults to None.

    Returns:
        (gpd.GeoDataFrame, gpd.GeoDataFrame): tuple of (edges, nodes) GeoDataFrames
    """

    # check if files exist
    flow_dir_fn = Path(flow_dir_fn)
    assert flow_dir_fn.exists()

    stations_fn = Path(stations_fn)
    assert stations_fn.exists()

    save_dir = Path(save_dir)
    if not save_dir.exists():
        print(f"passed save_dir does not exist, creating {save_dir}")
        save_dir.mkdir(parents=True)
    
    if elevation_fn is not None:
        elevation_fn = Path(elevation_fn)
        assert elevation_fn.exists()

    fdr = rxr.open_rasterio(flow_dir_fn, masked=True)
    band = fdr.sel(band=1)

    band_vicfmt = band

    reservoirs = gpd.read_file(stations_fn)
    reservoirs['geometry'] = gpd.points_from_xy(reservoirs['lon'], reservoirs['lat'])
    reservoirs.set_crs('epsg:4326', inplace=True)

    reservoir_location_raster = xr.full_like(band_vicfmt, np.nan)
    for resid, row in reservoirs.iterrows():
        reslat = float(row.lat)
        reslon = float(row.lon)

        rast_lat = reservoir_location_raster.indexes['y'].get_indexer([reslat], method="nearest")[0]
        rast_lon = reservoir_location_raster.indexes['x'].get_indexer([reslon], method="nearest")[0]

        reservoir_location_raster[rast_lat, rast_lon] = resid

    # convert all points to nodes. Use index value to identify
    G = gnx.GeoDiGraph()
    G.add_nodes_from(reservoirs.index)

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

    for node in G.nodes:
        resdata = reservoirs[reservoirs.index==node]
        
        x = float(resdata['lon'].values[0])
        y = float(resdata['lat'].values[0])

        idxx = band_vicfmt.indexes['x'].get_indexer([x], method="nearest")[0]
        idxy = band_vicfmt.indexes['y'].get_indexer([y], method="nearest")[0]
        
        # travel downstream until another node, np.nan or out-of-bounds is found, or if travelling in a loop

        visited = [(idxx, idxy)]
        current_pix = band_vicfmt.isel(x=idxx, y=idxy)

        attrs_n = {
            node: {
                'x': reservoirs['geometry'][node].x,
                'y': reservoirs['geometry'][node].y,
                'name': reservoirs['name'][node]
            }
        }
        nx.set_node_attributes(G, attrs_n)

        if not np.isnan(current_pix):
            END = False
            while not END:
                op = operations[int(current_pix)]
                new_idxy, new_idxx = np.array((idxy, idxx)) + np.array(op)
                idxy, idxx = new_idxy, new_idxx
                
                if (new_idxx, new_idxy) in visited:
                    # In a loop, exit
                    END=True
                    break
                else:
                    visited.append((new_idxx, new_idxy))

                current_pix = band_vicfmt.isel(x=new_idxx, y=new_idxy)
                if np.isnan(current_pix):
                    # NaN value found, exit loop
                    END=True
                    break

                try:
                    any_reservoir = reservoir_location_raster.isel(x=new_idxx, y=new_idxy)
                    if not np.isnan(any_reservoir):
                        # another reservoir found
                        G.add_edge(node, int(any_reservoir))
                        if dist_proj:
                            attrs_e = {
                                (node, int(any_reservoir)): {
                                    'length': reservoirs.to_crs(dist_proj)['geometry'][node].distance(reservoirs.to_crs(dist_proj)['geometry'][int(any_reservoir)])
                                }
                            }
                            nx.set_edge_attributes(G, attrs_e)
                        END = True
                        break
                except IndexError:
                    # Reached end
                    END=True

    G_gdf = gpd.GeoDataFrame(gnx.graph_edges_to_gdf(G))
    G_gdf_pts = gpd.GeoDataFrame(gnx.graph_nodes_to_gdf(G))

    # add elevation data to nodes
    if elevation_fn:
        elev = rxr.open_rasterio(elevation_fn, chunks='auto')

        G_gdf_pts['elevation'] = G_gdf_pts[['x', 'y']].apply(lambda row: float(elev.sel(x=row.x, y=row.y, method='nearest')), axis=1)
        G_gdf_pts.head()

    pts_save_fn = Path(save_dir) / 'rivreg_network_pts.shp'
    edges_save_fn = Path(save_dir) / 'rivreg_network.shp'
    
    G_gdf_pts.to_file(pts_save_fn)
    G_gdf.to_file(edges_save_fn)

    return G


def main():
    dist_proj = "+proj=eqdc +lon_0=-103.7988281 +lat_1=35.7127609 +lat_2=43.8942567 +lat_0=39.8035088 +datum=WGS84 +units=m +no_defs"
    flow_dir_fn = Path("data/colorado/basins/gunnison/ro/pars/fl.tif")
    stations_fn = Path("data/gunnison_reservoirs/gunnison_reservoirs_fdr_corrected_usgs_gages.csv")
    elevation_fn = Path("global_data/global_elevation_data/World_e-Atlas-UCSD_SRTM30-plus_v8.tif")
    save_fn = Path("data/gunnison_rivreg/rivreg_network.shp")
    save_fn.parent.mkdir(parents=True, exist_ok=True) # make parent directory if it doesn't exist

    assert flow_dir_fn.exists()
    assert stations_fn.exists()
    assert elevation_fn.exists()

    network, network_pts = generate_network(flow_dir_fn, stations_fn, dist_proj, elevation_fn)

    network.to_file(save_fn)
    network_pts.to_file(save_fn.with_name('rivreg_network_pts.shp'))


if __name__ == '__main__':
    main()
