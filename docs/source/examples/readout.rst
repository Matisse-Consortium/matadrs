.. role:: python(code)
   :language: python

=======
Readout
=======

This is a short guide how to use the readout functionality of the :ref:`matadrs` package.
First import :python:`pathlib` and then the :python:`ReadoutFits` from :python:`matadrs`.

.. code-block:: python

  from pathlib import Path

  from matadrs.utils.readout import ReadoutFits

Then specify a path to a (.fits)-file.

.. code-block:: python

  fits_file = Path("<fits_file>")

And finally we initialse the :python:`ReadoutFits` with the fits-file.

.. code-block:: python

  readout = ReadoutFits(file)

Now the (.fits)-file is loaded into the :python:`ReadoutFits`-class and the
values contained within the cards of the file can be accessed via the class's
attributes as tables.

.. code-block:: python

  print(readout.oi_wl)
  print(readout.oi_vis)
  print(readout.oi_vis.columns)
