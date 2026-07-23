# vesuvius_dem.tif — provenance and license

`vesuvius_dem.tif` is a Digital Elevation Model of Mt. Vesuvius (Italy) bundled
with UCMuon as the default sample terrain for Engine 6 (UCMuon Terrain) and the
MURAVES example in this directory.

| Property | Value |
|---|---|
| Product | SRTM GL1 (Shuttle Radar Topography Mission, Global 1 arc-second, ~30 m) |
| Bounding box | West 14.35 °E, South 40.76 °N, East 14.52 °E, North 40.90 °N |
| Grid | 612 x 504 pixels, EPSG:4326 (WGS84 lat/lon) |
| Elevation range | -8 m to 1262 m a.s.l. (Gran Cono summit) |
| Downloaded | 2026-07-20, via the OpenTopography Global DEM API |

## License

SRTM data are produced by NASA and USGS and are in the **public domain**
(U.S. Government work, no copyright). Redistribution of this file with UCMuon
is permitted without restriction.

## Citation

If you use this DEM in a publication, please cite the dataset:

> NASA Shuttle Radar Topography Mission (SRTM) (2013). *Shuttle Radar
> Topography Mission (SRTM) Global.* Distributed by OpenTopography.
> https://doi.org/10.5069/G9445JDF

and acknowledge OpenTopography (their requested text):

> This material is based on data services provided by the OpenTopography
> Facility with support from the National Science Foundation under NSF Award
> Numbers 1948997, 1948994 and 1948857.

## Re-downloading / other sites

The same file can be regenerated from the GUI (Terrain tab, Section 2,
"Auto-download (OpenTopography)") with the bounding box above, or from the
command line with a free personal API key from
<https://opentopography.org/requestapi>:

```bash
curl -o vesuvius_dem.tif "https://portal.opentopography.org/API/globaldem?demtype=SRTMGL1&south=40.76&north=40.90&west=14.35&east=14.52&outputFormat=GTiff&API_Key=YOUR_API_KEY"
```

Note: COP30 (Copernicus GLO-30) is an alternative 30 m product with slightly
better quality on steep slopes, but its license requires shipping the ESA/
Airbus license text alongside any redistributed copy, so the public-domain
SRTM GL1 is used for the bundled sample.
