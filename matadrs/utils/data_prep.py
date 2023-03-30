from pathlib import Path
from typing import Optional, List

import numpy as np
from astropy.table import Table, vstack
from astropy.table.np_utils import TableMergeError

from .readout import ReadoutFits

__all__ = ["DataPrep"]


class DataPrep:
    """A class that unifies the ReadoutFits class's functionality for
    multiple (.fits)-files

    It sources the data from the individual
    (.fits)-files and merges them into unified tables

    Parameters
    ----------
    fits_files: List[Path]
        The (.fits)-files from which data is sourced
    flux_files: List[Path], optional
        Additional flux files, that replace the flux data for a files
        that might be missing it

    Attributes
    ----------
    nfits: int
        The number of the read in (.fits)-files
    oi_wl
    oi_array
    oi_flux
    oi_t3
    oi_vis
    oi_vis2

    Methods
    -------
    """

    def __init__(self, fits_files: List[Path],
                 flux_files: Optional[List[Path]] = None) -> None:
        """Initialises the class"""
        self.fits_files = fits_files

        if flux_files is None:
            self.flux_files = [None]*len(fits_files)
        else:
            self.flux_files = flux_files

        self.nfits = len(self.fits_files)

        self._oi_wl, self._oi_array = [None]*2
        self._oi_flux, self._oi_t3 = [None]*2
        self._oi_vis, self._oi_vis2 = [None]*2

        self.readouts = list(map(ReadoutFits, self.fits_files))

    def __repr__(self):
        """The DataHandler class' representation"""
        return "DataHandler contains information on the following (.fits)-files\n"\
            + f"{'':-^50}\n"\
            + '\n'.join([fits_file.name for fits_file in self.fits_files])

    def __str__(self):
        """The DataHandler class' string printout"""
        return self.__repr__()

    @property
    def longest_entry(self):
        """The longest entry of all the rows. Fetched from the 'oi_wl'-tables"""
        return np.max([readout.longest_entry for readout in self.readouts])

    @property
    def oi_wl(self):
        """Gets the unified wavelength table"""
        if self._oi_wl is None:
            self._oi_wl = self._set_table_attribute("oi_wl")
        return self._oi_wl

    @property
    def oi_array(self):
        """Gets the unified wavelength table"""
        if self._oi_array is None:
            self._oi_array = self._set_table_attribute("oi_array")
        return self._oi_array

    @property
    def oi_flux(self):
        """Fetches the flux table"""
        # TODO: Check how the flux can be unified and later how to handle empty/None
        # rows!!!
        if self._oi_flux is None:
            self._oi_flux = self._set_table_attribute("oi_flux")
        return self._oi_flux

    @property
    def oi_vis(self):
        """Fetches the visibility table"""
        if self._oi_vis is None:
            self._oi_vis = self._set_table_attribute("oi_vis")
        return self._oi_vis

    @property
    def oi_vis2(self):
        """Fetches the squared visibility table"""
        if self._oi_vis2 is None:
            self._oi_vis2 = self._set_table_attribute("oi_vis2")
        return self._oi_vis2

    @property
    def oi_t3(self):
        """Fetches the closure phase table"""
        if self._oi_t3 is None:
            self._oi_t3 = self._set_table_attribute("oi_t3")
        return self._oi_t3

    def _set_table_attribute(self, attribute_name: str):
        """Sets and unifies the input tables for the attributes given"""
        if len(self.readouts) > 1:
            try:
                return vstack([getattr(readout, attribute_name)
                               for readout in self.readouts])
            except TableMergeError:
                # TODO: implement padding here
                ...
        else:
            return getattr(self.readouts[0], attribute_name)

    # TODO: Check how to pad the arrays so they all have length 121
    def pad_column(self, table: Table, columns: List[str]):
        """Pads the individual columns with 'np.nan' up to a certain length"""
        for colname in table.columns:
            for row in table[colname]:
                ...
