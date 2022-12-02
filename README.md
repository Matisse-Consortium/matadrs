# Matadrs
The (modular) MATISSE Automated Data Reduction Software

## Features
### Raw-data fetching from ESO-Archive
> _To be implemented_

### Data-reduction for MATISSE
The script `matadrs.py` gives easy access of functionality from the following scripts
(found under `libs/reduction/`)
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

### Mat-tools (required)
Will be automatically installed, for in-depth information see
[mat-tools](https://gitlab.oca.eu/MATISSE/tools/-/tree/master/mat_tools)

## Matadrs
To install this package (locally) run
```
pip3 install -e .
```

### Virtual Enviroment Setup (optional)
There are multiple ways to create a virtual enviroment in python, such as
[pyenv](https://towardsdatascience.com/managing-virtual-environment-with-pyenv-ae6f3fb835f8),
which is an extension of the python (built-in) package
[virtualenv](https://learnpython.com/blog/how-to-use-virtualenv-python/) or other tools
such as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html).

## Usage

## Ressources

## Contributors
[M. B. Scheuck](https://github.com/MBSck)
