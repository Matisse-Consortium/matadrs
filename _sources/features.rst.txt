.. role:: bash(code)
   :language: bash

========
Features
========

For working examples of all the features, and more explanations on the usage
see :ref:`Getting Started` and for an in depth view of the individual code,
see :ref:`api`.

Raw Data Retrieval
==================

Automatic download of the data for one or more specific night(s) from the ESO-archives.

*Yet to be implemented*

Reduction
=========

The reduction can be done with simply using the function
:func:`matadrs.matadrs.reduction_pipeline`.

The function :func:`matadrs.matadrs.reduction_pipeline` combines functionality
from the following modules (found under the :ref:`matadrs.reduction <reduction_pkg>`
subpackage).

* :mod:`matadrs.reduction.reduce`
* :mod:`matadrs.reduction.calibrate`
* :mod:`matadrs.reduction.average`
* :mod:`matadrs.reduction.merge`

These modules take care of the reduction, calibration of the fluxes and visibilities,
averaging all datasets together, and finally merging all of them into one final,
oifits-compliant (.fits)-file.

Readout
=======

This enables the user to read in a (.fits)-file and get the data in the form of a
`astropy.table <https://docs.astropy.org/en/stable/table/index.html>`_ for every
header in the (.fits)-file (:bash:`oi_vis`, :bash:`oi_t3`, etc.),
as long as the file is in the `oifits <https://oifits.org/>`_ format.

Plotting
========

There is additional plotting functionality included in :bash:`matadrs` in order to
immediately visualise the products of the data-reduction pipeline at each individual
step (see :ref:`reduction`).
