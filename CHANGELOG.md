# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- SMAP L4 upgraded to v8.

### Removed
- Support for OpenDAP access to SMAP.

### Fixed

- Added initialization to exposure forecast file names to avoid overwriting.

## [2.1] - 2024-9-3

### Added

- Slope threshold parameter to CLI.

### Changed

- Reduced size of exposure table.
- Streamlined conda environment for global landslide nowcast. The PFDF model will use a separate environment.
- Increased default cache size to last 7 days of input data.
- Switched to daily GIS data on IMERG server.
- Recalculated GEOS 99th percentile from recent (2022-2023) forecast data.
- Changed seismic preconditioning variable to PGA.
- Trained new model with additional landslide inventories.

### Removed
- Support for OpenDAP access to IMERG.
- Snow mass variable, as XGBoost did not use it in v2.1

### Fixed

- Switched default IMERG version to 07B.

## [2.0] - 2023-6-29

### Changed

- Replaced land mask with a file based on the MOD44W global water mask. This is combined with the existing mask from SMAP L4. To use the new mask file, users should download [static.zip](https://gpm.nasa.gov/sites/default/files/data/landslides/static.zip) again.
- Extended SMAP file search to 5 days to deal with cases where latency exceeds the typical 3 days.
- SMAP L4 upgraded to v7.
- Reduced memory usage by exposure function.

### Fixed

- Replaced the post-fire debris flow model, which had stopped working due to changes to the Google Earth Engine API. (#8).
- Forecasted exposure.

## [2.0-beta] - 2022-7-5

### Added

- Option to set number of threads for predictions.
- Option to download data instead of using OpenDAP.
- Updated citations and metadata.

### Fixed

- Unzipping during installation.
- Location of subsetted geotiffs.

## [2.0-alpha] - 2022-4-23

### Added

- Initial release of LHASA 2.0.

## [1.1.1] - 2021-08-24

### Changed

- Updated ARI threshold.
- Minor changes for ease of use.

## [1.1.0] - 2020-4-2

### Changed

- Updated landslide susceptibility map.
- Switched from total precipitation to rainfall as an input.

## [1.0.0] - 2019-3-13

### Added

- Initial publication of R code for NRT landslide hazard mapping.
