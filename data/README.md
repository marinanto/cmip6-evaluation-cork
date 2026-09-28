# Input data

This directory contains the climate datasets and station observations used in the County Cork climate analysis.

## Gridded climate data

The `new_eraland/` directory contains the ERA5-Land reference dataset and CMIP6 GCM simulations used in the analysis.

| File | Dataset | Period | Use in the analysis |
|---|---|---|---|
| `era5_land_cork_1985_2025_final.nc` | ERA5-Land | 1985-2025 | Reference dataset for GCM evaluation, bias correction and validation |
| `cnrm_ref_final.nc` | CNRM | 1985-2014 | Baseline GCM evaluation |
| `hadgem_ref_final.nc` | HadGEM | 1985-2014 | Baseline GCM evaluation |
| `ecearth_ref_final.nc` | EC-Earth | 1985-2014 | Baseline GCM evaluation |
| `cnrm45_final.nc` | CNRM, SSP2-4.5 | 2015-2025 | Future validation |
| `cnrm85_final.nc` | CNRM, SSP5-8.5 | 2015-2025 | Future validation |

### ERA5-Land

ERA5-Land data were obtained from the Copernicus Climate Change Service (C3S) Climate Data Store. ERA5-Land provides gridded land-surface information at approximately 9 km native resolution and includes near-surface air temperature.

**Dataset:** ERA5-Land hourly data from 1950 to present  
**DOI:** 10.24381/cds.e2161bac  
**Licence:** CC BY  
**Provider:** Copernicus Climate Change Service (C3S), implemented by ECMWF

Source:

https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land

The ERA5-Land dataset is described by Muñoz-Sabater et al. (2021):

> Muñoz-Sabater, J., Dutra, E., Agustí-Panareda, A., Albergel, C., Arduini, G., Balsamo, G., Boussetta, S., Choulga, M., Harrigan, S., Hersbach, H., Martens, B., Miralles, D. G., Piles, M., Rodríguez-Fernández, N. J., Zsoter, E., Buontempo, C., & Thépaut, J.-N. (2021). ERA5-Land: A state-of-the-art global reanalysis dataset for land applications. *Earth System Science Data*, 13, 4349-4383. https://doi.org/10.5194/essd-13-4349-2021

### CMIP6 climate projections

The CNRM, HadGEM and EC-Earth simulations were obtained through the Copernicus Climate Change Service (C3S) Climate Data Store CMIP6 climate projections dataset.

**Dataset:** CMIP6 climate projections  
**DOI:** 10.24381/cds.c866074c  
**Provider:** Copernicus Climate Change Service (C3S), Climate Data Store  
**Experiments used:** historical/reference simulations and SSP2-4.5 and SSP5-8.5 scenarios

Source:

https://cds.climate.copernicus.eu/datasets/projections-cmip6

The CMIP6 dataset provides daily and monthly climate projections from multiple CMIP6 experiments and GCMs and is distributed through the C3S Climate Data Store.

## Station observations

The `stations/` directory contains monthly mean temperature observations from four Met Éireann stations in County Cork:

| File | Station |
|---|---|
| `mly575.xlsx` | Moore Park |
| `mly775.xlsx` | Sherkin Island |
| `mly1075.xlsx` | Roches Point |
| `mly3904.xlsx` | Cork Airport |

The station data are used as an independent observational dataset for validation of the raw and bias-corrected CNRM simulations during 2015-2025.

The station files contain the following variables:

- `year` - observation year
- `month` - observation month
- `meant` - monthly mean temperature (°C)

### Station data source

The station observations were obtained from Met Éireann, the Irish Meteorological Service. Met Éireann provides historical climate observations at hourly, daily and monthly temporal resolutions, including temperature observations from its national observing network.

Source:

https://www.met.ie/climate/available-data

Monthly station data:

https://www.met.ie/climate/available-data/monthly-data

Met Éireann National Observing Network:

https://www.met.ie/climate/the-national-observing-network

## Data use and attribution

The datasets in this directory are retained with the repository to support reproducibility of the analysis. The original data providers retain the applicable intellectual property rights and dataset-specific terms of use.

Users reusing these datasets should cite the original data sources and comply with their applicable licences and attribution requirements.

For ERA5-Land, the applicable C3S dataset licence is CC BY. For CMIP6 data, users should follow the CMIP6 data access terms specified by the Copernicus Climate Change Service. Met Éireann data should be used and attributed according to the terms specified by Met Éireann.

## Reproducibility

The filenames and directory structure correspond to the inputs expected by:

`R/cork_analysis_9km.R`