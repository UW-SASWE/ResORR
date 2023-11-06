# ResORR - Reservoir Operations driven River Regulation

Model for estimating river regulation due to reservoir operations using the [Reservoir Assessment Tool 3.0](https://github.com/UW-SASWE/RAT/).

## `rrview` - interactive dashboard for exploring ResORR data

Run the dashboard using the following command - 
```bash
panel serve src/rrview/app.py --show --args <river-regulation-data-file> -npts <path-to-network-points-file> -nlnk <path-to-network-links-file>
```