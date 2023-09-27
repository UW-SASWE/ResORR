# River Regulation

Model for estimating river regulation due to dams using the [Reservoir Assessment Tool 3.0](https://github.com/UW-SASWE/RAT/).

## `rrview` - interactive dashboard for exploring river regulation data

Run the dashboard using the following command - 
```bash
panel serve src/rrview/app.py --show --args <river-regulation-data-file> -npts <path-to-network-points-file> -nlnk <path-to-network-links-file>
```

Example: 
```bash
panel serve src/rrview/app.py --show --args /tiger1/pdas47/work/2023_01_24-river-regulation/data-era5_2010_2021/regulation/regulation_data.insitu.regulated.ERA5-2010_2021.nc -npts /tiger1/pdas47/work/2023_01_24-river-regulation/data-cumberland/cumberland_rivreg/cumberland_rivreg_pts.geojson -nlnk /tiger1/pdas47/work/2023_01_24-river-regulation/data-cumberland/cumberland_rivreg/cumberland_rivreg.geojson
```