from pathlib import Path
from typing import Optional, List

import astropy.units as u
import numpy as np
import pkg_resources
from astropy.table import Table, vstack
from astropy.table.np_utils import TableMergeError

from .readout import ReadoutFits

__all__ = []


DATA_DIR = Path(pkg_resources.resource_filename("matadrs", "data"))


# TODO: Keep all in Table format as long as possible but switch to np.ndarrays .when needed
# TODO: Find way to unify tables so that all sub-tables are still differentiable and
# recognisable, but also that the unfied values are easily accesible
# TODO: Implement the np.nan functionality in plotting -> For fits plotting and plotter
# class. Implement fits plotting in plotter class as well
# TODO: For model fitting. Save the infos as astropy.Table
# TODO: Implement all the properties from ReadoutFits in this class as well
# TODO: Maybe rename the class to better fit its purpose
# TODO: Make a submodule that contains the plot.py, readout.py and data_prep.py
# TODO: Implement the flux-files somehow
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

        self.readouts = list(map(ReadoutFits, self.fits_files))

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
            self._oi_wl = self._set_table_attribute("oi_wl")
        return self._oi_wl

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
        """Sets the Table properties"""
        if len(self.readouts) > 1:
            try:
                return vstack([getattr(readout, attribute_name)\
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

    def get_data_for_wavelength(self, table: Table, wavelengths: u.um) -> Table:
        """Fetches the corresponding data for the sought wavelength"""
        ...

    # def get_data_for_wavelength(self, data: Union[Quantity, np.ndarray],
                                # wl_poly_indices: List) -> List:
        # """Fetches data for one or more wavelengths from the nested arrays. Gets the
        # corresponding values by index from the nested arrays (baselines/triangle)

        # Parameters
        # ----------
        # data: astropy.units.Quantity | numpy.ndarray
            # The data for every baseline/triangle
        # wl_poly_indices: List
            # The polychromatic indices of the wavelength solution. This has to be a doubly
            # nested list

        # Returns
        # --------
        # data4wl: List
        # """
        # # NOTE: Right now the data is immediately averaged after getting taken. Maybe
        # # change this for the future
        # polychromatic_data_averaged = []
        # for dataset in data:
            # data4wl = []
            # for wl_indices in wl_poly_indices:
                # data4wl_poly_index = []
                # for wl_index in wl_indices:
                    # array_wl_slice = u.Quantity([array[wl_index] for array in dataset])
                    # data4wl_poly_index.append(array_wl_slice)
                # data4wl.append(u.Quantity(data4wl_poly_index))
            # averaged_dataset_slice = self.average_polychromatic_data(data4wl)
            # polychromatic_data_averaged.append(averaged_dataset_slice)
        # return [u.Quantity(dataset4wl) for dataset4wl in polychromatic_data_averaged]

    # def average_polychromatic_data(self, polychromatic_data: Quantity):
        # """Fetches and then averages over polychromatic data. Iterates over the
        # polychromatic wavelength slices and then takes the mean of them

        # Parameters
        # ----------
        # polychromatic_data: astropy.units.Quantity
            # The polychromatic data slices of wavelengths in one window
        # """
        # return u.Quantity([np.mean(data_slice, axis=0)\
                           # for data_slice in polychromatic_data])


if __name__ == "__main__":
    fits_files = ["HD_163296_2019-03-23T08_41_19_N_TARGET_FINALCAL_INT.fits"]
                  # "HD_163296_2019-03-23T08_41_19_L_TARGET_FINALCAL_INT.fits",
                  # "HD_163296_2019-05-06T08_19_51_L_TARGET_FINALCAL_INT.fits"]
    fits_files = [DATA_DIR / "tests" / fits_file for fits_file in fits_files]
    data_prep = DataPrep(fits_files)
    breakpoint()
