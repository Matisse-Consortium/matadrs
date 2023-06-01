.. role:: python(code)
   :language: python

========
Plotting
========

This is a short guide how to use the plotting functionality of the :ref:`matadrs` package.
First import :python:`pathlib` and then the :python:`Plotter` from :python:`matadrs`.

.. code-block:: python

  from pathlib import Path

  from matadrs.utils.plot import Plotter

Then specify a path to one or more valid (.fits)-file(s).

.. code-block:: python

  fits_file = Path("<fits_file>")

And finally we initialse the :python:`Plotter` with the fits-file(s).

.. code-block:: python

  plotter = Plotter(fits_file)

Now one can choose from multiple options for the plots, some of which are presented
in the following (for more options see the api's documentation).

.. code-block:: python

  plotter.add_uv().add_vis(corr_flux=True).add_cphases(unwrap=True)
  plotter.plot(save=True, format="pdf")

As plotter returns :python:`self` with every function, one can call all the function including
:python:`plot` in the same line.

.. code-block:: python

  plotter.add_uv().plot(save=False)
