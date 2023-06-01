.. role:: python(code)
   :language: python

.. role:: bash(code)
   :language: bash

======================
Executing the Pipeline
======================

This is a short guide how to use the :bash:`matadrs` data-reduction pipeline.
First import :python:`pathlib` and the :func:`matadrs.matadrs.reduction_pipeline`
from :bash:`matadrs`.

.. code-block:: python

  from pathlib import Path

  from matadrs import reduction_pipeline

Now the paths to the directories with the ESO-archive data, the raw-directories
, :python:`raw_dirs`, need to be specified. This can be multiple directories,
then the pipeline will reduce all these data. For instance, a data directory
and then multiple nights. Also the output directories, :python:`product_dirs`,
need to be given as well.

.. code-block:: python

  # Specify the path to the directory containing the data
  data_dir = Path("<data_dir>")

  # Specify the raw-directory, containing the raw data
  observation_dirs = ["<night1>", "<night2>"]
  raw_dirs = list(map(lambda x: data_dir / "raw" / x,
                      observation_dirs))

  # Specify the product-directory, to contain the product data/that contains reduced,
  # calibrated or averaged data, to be further processed
  product_dirs = list(map(lambda x: Path(str(x).replace("raw", "product")),
                          raw_dirs))

Then with all that give, the :ref:`data-reduction pipeline <matadrs.matadrs.reduction_pipeline>`
can be called

.. code-block:: python

  reduction_pipeline(raw_dirs, product_dirs,
                     overwrite=True, do_reduce=True,
                     do_calibrate=True, do_average=True,
                     do_merge=True)

All the settings
