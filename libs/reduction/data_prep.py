from pathlib import Path
from typing import Optional, List

import numpy as np
import astropy.units as u
from astropy.table import Table, vstack
from astropy.table.np_utils import TableMergeError

from readout import ReadoutFits


DATA_DIR = Path(__file__).parent.parent.parent / "data"


# TODO: Keep all in Table format as long as possible but switch to np.ndarrays when needed
# TODO: Find way to unify tables so that all sub-tables are still differentiable and
# recognisable, but also that the unfied values are easily accesible
# TODO: Implement the np.nan functionality in plotting -> For fits plotting and plotter
# class. Implement fits plotting in plotter class as well
# TODO: For model fitting. Save the infos as astropy.Table
# TODO: Implement all the properties from ReadoutFits in this class as well
# TODO: Maybe rename the class to better fit its purpose
class DataPrep:
    """"""

    def __init__(self, fits_files: List[Path],
                 flux_files: Optional[List[Path]] = None) -> None:
        """Initialises the class"""
        self.fits_files = fits_files
        if flux_files is None:
            self.flux_files = [None]*len(fits_files)
        else:
            self.flux_files = flux_files

        self._oi_wl = None
        self._oi_flux, self._oi_t3 = None, None
        self._oi_vis, self._oi_vis2 = None, None

        self.readouts = [ReadoutFits(fits_file) for fits_file in self.fits_files]

        # TODO: Implement this in plotter
        # self.uv_coords = self._merge_simple_data("uvcoords")
        # self.uv_coords_cphase = self._merge_simple_data("uvcoords_cphase")
        # self.baselines_cphase = self._merge_simple_data("baselines_cphase")
        # self.telescope_info = self._merge_simple_data("telescope")
        # self.baselines = self._merge_simple_data("baselines")

    def __repr__(self):
        """The DataHandler class' representation"""
        return "DataHandler contains information on the following (.fits)-files:\n"\
                + '\n'.join([fits_file.name for fits_file in self.fits_files])

    def __str__(self):
        """The DataHandler class' string representation"""
        return self.__repr__()

    @property
    def longest_entry(self):
        """The longest entry of all the rows. Fetched from the 'oi_wl'-tables"""
        return np.max([readout.longest_entry for readout in self.readouts])

    @property
    def oi_wl(self):
        """Gets the unified wavelength table"""
        if self._oi_wl is None:
            try:
                self._oi_wl = vstack([readout.oi_wl for readout in self.readouts])
            except TableMergeError:
                # TODO: implement padding here
                ...
        return self._oi_wl

    @property
    def oi_flux(self):
        """Fetches the flux table"""
        # TODO: Check how the flux can be unified and later how to handle empty/None
        # rows!!!
        if self._oi_flux is None:
            try:
                self._oi_wl = vstack([readout.oi_wl for readout in self.readouts])
            except TableMergeError:
                # TODO: implement padding here
                ...
        return self._oi_flux

    @property
    def oi_vis(self):
        """Fetches the visibility table"""
        if self._oi_vis is None:
            try:
                self._oi_vis = vstack([readout.oi_vis for readout in self.readouts])
            except TableMergeError:
                # TODO: implement padding here
                ...
        return self._oi_vis

    @property
    def oi_vis2(self):
        """Fetches the squared visibility table"""
        if self._oi_vis2 is None:
            try:
                self._oi_vis2 = vstack([readout.oi_vis2 for readout in self.readouts])
            except TableMergeError:
                # TODO: implement padding here
                ...
        return self._oi_vis2

    @property
    def oi_t3(self):
        """Fetches the closure phase table"""
        if self._oi_t3 is None:
            try:
                self._oi_t3 = vstack([readout.oi_t3 for readout in self.readouts])
            except TableMergeError:
                # TODO: implement padding here
                ...
        return self._oi_t3

    def get_data_for_wavelength(self, table: Table, wavelengths: u.um) -> Table:
        """Fetches the corresponding data for the sought wavelength"""
        ...

    # TODO: Check how to pad the arrays so they all have length 121
    def pad_column(self, table: Table, columns: List[str]):
        """Pads the individual columns with 'np.nan' up to a certain length"""
        for colname in table.columns:
            for row in table[colname]:
                ...


if __name__ == "__main__":
    fits_files = ["HD_163296_2019-03-23T08_41_19_N_TARGET_FINALCAL_INT.fits",
                  "HD_163296_2019-03-23T08_41_19_L_TARGET_FINALCAL_INT.fits"]
                  # "HD_163296_2019-05-06T08_19_51_L_TARGET_FINALCAL_INT.fits"]
    fits_files = [DATA_DIR / "tests" / fits_file for fits_file in fits_files]
    data_prep = DataPrep(fits_files)
    breakpoint()
