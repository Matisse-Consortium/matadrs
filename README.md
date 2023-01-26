# Matadrs
The (modular) MATISSE Automated Data Reduction Software

## Features
### Raw-data fetching from ESO-Archive
> _To be implemented_

### Data-reduction for MATISSE
The script `matadrs.py` gives easy access of functionality from the following scripts
(found under `matadrs/reduction/`)
* `reduce.py`
* `calibrate.py`
* `average.py`
* `merge.py`

## Installation
It is recommended to keep your Python3 installation clean by installing packages to a
[virtual-environment](#virtual-enviroment-setup)

### MATISSE-pipeline (required)
First the MATISSE-pipeline needs to be installed.<br>
For this follow the instructions found here:

## Matadrs
To install this package (locally, in development mode) run
```
pip3 install -e .
```
For a permanent installation run
```
pip3 install matadrs git+ssh://git@github.com:MBSck/matadrs.git
```

## Usage

## Ressources

## Credit
* J. Varga for his utmost helpful input and the `fluxcal` and `avg_oifits` scripts
* R. van Boekel, J. Varga, A. Matter for their calibrator and object catalogues as well as
  the [JSDC-catalog ](https://cdsarc.cds.unistra.fr/viz-bin/cat/II/300)
* F. Millour for the `calib_BCD2` script and J. Varga for his rework of the same
* M. Neeser for the `getdata_ESO_archive` script
