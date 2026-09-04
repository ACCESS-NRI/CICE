# JRA55-do forcing generation

`generate_jra55do_forcing.py` regrids one year of JRA55-do (input4MIPs v1.6.0) atmospheric
forcing onto a CICE standalone grid, producing a single annual file in the format CICE's
`ice_forcing.F90` `JRA55(do)` option expects (see **JRA55_files**/**JRA55_data** in
`cicecore/cicedyn/general/ice_forcing.F90`, and `dg_forcing.rst` in the Developer Guide).

Unlike `../jra55_datasets/interp_jra55_ncdf_bilinear.py` (which targets the older, raw JRA55
reanalysis file layout), this script reads directly from the JRA55-do input4MIPs source tree
and only requires the destination grid's MOM6 supergrid file.

## Setup

```
conda env create -f environment.yml
conda activate jra55do_datasets
```

## Usage

```
python generate_jra55do_forcing.py --year=<year> --hgrid-filename=<path-to-supergrid-file> \
    --output-filename=<path-to-output-file>
```

By default this reads from the JRA55-do v1.6.0 source tree at
`/g/data/qv56/replicas/input4MIPs/CMIP6Plus/OMIP/MRI/MRI-JRA55-do-1-6-0` (NCI Gadi). Pass
`--jra55do-dir` to point at a different copy of the source data.

Run `python generate_jra55do_forcing.py -h` for the full list of options.

## Method

Seven fields are read from the raw annual JRA55-do source files, regridded (bilinear) onto
the destination grid's T-cell centres, renamed to CICE's expected variable names, and
combined into one file:

| JRA55-do variable | CICE variable | Frequency                              |
|---|---|---|
| `tas`   | `airtmp` | instantaneous, 8x daily |
| `uas`   | `wndewd` | instantaneous, 8x daily |
| `vas`   | `wndnwd` | instantaneous, 8x daily |
| `huss`  | `spchmd` | instantaneous, 8x daily |
| `rsds`  | `glbrad` | 3-hour average |
| `rlds`  | `dlwsfc` | 3-hour average |
| `prsn`  | `ttlpcp` | 3-hour average |

No intermediate Repeat-Year-Forcing (RYF) step is needed for a single real calendar year -
the raw annual source files already have the correct record count and order.

**Precipitation note**: `ttlpcp` is built from JRA55-do's `prsn` (snowfall) alone, not
`prsn + prra` (total precipitation). CICE's JRA55(do) reader has no rain/snow partitioning
and adds `ttlpcp` directly to `fsnow` (see `ice_forcing.F90:JRA55_data`), so including
rainfall would incorrectly convert it to snow mass.

## Output

A single NetCDF file with 3-hourly (8x daily) records for the requested year, on the
destination grid's tracer (T-cell) grid. The generating command, script git commit (if
applicable), and input file paths are recorded in the output file's global attributes.
