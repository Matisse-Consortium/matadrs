# Matadrs
The (modular) MATISSE Automated Data Reduction Software (Matadrs) is a tool written to
encompass all the individual steps and files needed to procure and reduce data generally
and easily in as few steps as possible with the least knowledge required. To make it
possible for seasoned MATISSE users to quickly and easily reduce data and giving others
who are interested in MATISSE data the opportunity to do so as well.

## Features
### Raw-data fetching from ESO-Archive
> _To be implemented_

### Data-reduction for MATISSE
This is a wrapper around various other scripts and the MATISSE-data reduction pipeline
that can be used to reduce the raw-data from the ESO-Archive.
The script `matadrs.py` uses functionality from the following main-scripts (found under `matadrs/reduction/`)
* `reduce.py`
* `calibrate.py`
* `average.py`
* `merge.py`

That take care of the reduction, then calibrate the fluxes and visibilities and average
all the datasets together to finally merge them all into one final (.fits)-file.<br>
There is also additional functionality contained in the `matadrs/utils` directory, for
instance, the  `readout`- and `plot`-scripts.

## Installation
It is recommended to keep your Python3 installation clean by installing packages within a
[virtual-environment](#virtual-enviroment-setup)

### MATISSE-pipeline (required)
Follow the instructions found here: https://www.eso.org/sci/software/pipelines/

### Matadrs
To install this package (in development mode/editable mode) run
```
pip3 install -e .
```
For a permanent/non-editable installation run
```
pip3 install git+ssh://git@github.com:MBSck/matadrs.git
```

## Usage

For detailed usage see the `example`-directory.<br>
_Soon the to be upcoming docs as well!_

## Credit
* J. Varga for his utmost helpful input and his `fluxcal`- and `avg_oifits`-scripts
* R. van Boekel, J. Varga, A. Matter for their calibrator- and object catalogues as well as
  the [JSDC-catalog](https://cdsarc.cds.unistra.fr/viz-bin/cat/II/300)(JMMC Stellar Diameter Catalog) from Lafrasse et al.
* F. Millour for the `calib_BCD2`-script and J. Varga for his rework of the same
* M. Neeser for the `getdata_ESO_archive`-script
