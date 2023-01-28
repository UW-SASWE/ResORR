import ee
import geopandas as gpd
from pathlib import Path

# authorize using service key
credentials = ee.ServiceAccountCredentials('saswe-laura@pdas-river-regulation.iam.gserviceaccount.com', 'secrets/pdas47-river-regulation-ca5e5302cda7.json')


ee.Initialize(credentials)

def main():
    # download a single asset
    asset = ee.FeatureCollection("projects/pdas47-river-regulation/assets/rat3_reservoir_shp")
    asset = asset.getInfo()

    gdf = gpd.GeoDataFrame.from_features(asset)

    gdf.to_file('data/rat3_reservoir_shp.geojson', driver='GeoJSON')

if __name__ == '__main__':
    main()