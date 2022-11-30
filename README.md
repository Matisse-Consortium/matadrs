# Matadrs
The (modular) MATISSE Automated Data Reduction Software

## Table of Contents
* [Features](#features)
* [Installation](#installation)
* [Usage](#usage)
* [Ressources](#ressources)
* [Contributors/Contact](#contributors)

## Features
* `reduction.py`
* `calibration.py`
* `average.py`
* `merge.py`

## Installation
> It is recommended to keep your Python3 installation clean by installing packages to a
virtual-environment (for this see [Virtual Enviroment Setup](#virtual-enviroment-setup))

### MATISSE-pipeline (required)

### Mat-tools (optional)
This should be automatically installed with the `requirements.txt`,
however, if you want to install a newer version follow the instructions from
[installation-instructions](https://gitlab.oca.eu/MATISSE/tools/-/wikis/home)

> :warning: Sometimes there occurs an error while installing the `wxPython`-wheels,
required for the `mat-tools` package, this can be fixed by manually installing
the missing package (the error should tell you which it is)

## Matadrs
After installing the requirements (see [MATISSE-pipeline](#matisse-pipeline)
and (possibly) manually installing the [Mat-tools](#mat-tools) package)
you can install the remaining dependencies via executing,
```
pip install requirements.txt
```

## Usage

## Ressources
### Virtual Enviroment Setup
There are multiple ways to create a virtual enviroment in python, such as
[pyenv](https://towardsdatascience.com/managing-virtual-environment-with-pyenv-ae6f3fb835f8),
which is an extension of the python (built-in) package
[virtualenv](https://learnpython.com/blog/how-to-use-virtualenv-python/) or other tools
such as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html).

## Contributors
[M. B. Scheuck](https://github.com/MBSck)
