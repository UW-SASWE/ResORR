import ee

ee.Initialize()

from rat.core.sarea.sarea_cli_l8 import sarea_l8
import geopandas as gpd
from pathlib import Path

res_extent_fn = "/water2/pdas47/2023_01_24-river-regulation/data-cumberland/cumberland-stations/cumberland-reservoirs.geojson"
data_dir = Path("/water2/pdas47/2023_01_24-river-regulation/data-cumberland/dels")
start = "2018-01-01"
end = "2022-12-31"

res_extent = gpd.read_file(res_extent_fn)
l8_dir = data_dir / "l8"
l8_dir.mkdir(exist_ok=True, parents=True)

for i, row in res_extent.iterrows():
    name = row["name"]
    print(f"!!!! Processing {name}")

    geom = row.geometry
    sarea_l8(name, geom, start, end, l8_dir)