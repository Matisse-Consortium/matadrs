from typing import Dict, List

import numpy as np
import astropy.units as u
from astropy.io import fits
from pathlib import Path
from astropy.table import Table


DATA_DIR = Path(__file__).parent.parent.parent.parent.parent / "data"


class ReadoutFits:
    """All functionality to work with '.oifits/.fits'-files"""
    def __init__(self, fits_file: Path) -> None:
        self.fits_file = Path(fits_file)

        self._wl = None
        self._telescope_to_station_names = None

        with fits.open(self.fits_file) as hdul:
            self.primary_header = hdul[0].header

        self.target_name = self.primary_header["OBJECT"]
        self.ra, self.dec = self.primary_header["RA"], self.primary_header["DEC"]

    # def __repr__(self):
        # """The class's representation"""
        # # TODO: Make this better
        # return f"The instance is maintaining the following files:\n"\
               # f"{', '.join([fits_file.name for fits_file in self.fits_file])}"

    # def __str__(self):
        # """The class's print"""
        # return self.__repr__()

    @property
    def wl(self):
        """Fetches the wavelengths"""
        if self._wl is None:
            self._wl = self.get_table_for_fits("oi_wl")["EFF_WAVE"].to(u.um)
        return self._wl

    @property
    def telescopes_to_station_names(self) -> List[Dict]:
        """Makes a dict of the station indicies and their station names"""
        if self._telescope_to_station_names is None:
            oi_array = self.get_table_for_fits("oi_array")
            self._telescope_to_station_names = dict(zip(oi_array["STA_INDEX"],
                                                        oi_array["TEL_NAME"]))
        return self._telescope_to_station_names

    @property
    def bcd_configuration(self):
        """Gets the BCD-configuration"""
        return ["-".join([self.primary_header["HIERARCH ESO INS BCD1 ID"],
                          self.primary_header["HIERARCH ESO INS BCD2 ID"]])]

    # def get_split_uvcoords(self) -> np.ndarray:
        # """Splits a 2D-np.array into its 1D-components and returns the u- and
        # v-coords seperatly"""
        # uvcoords = self.get_uvcoords()
        # return np.array([item[0] for item in uvcoords]),\
                # np.array([item[1] for item in uvcoords])

    def get_table_for_fits(self, header: str) -> List[Table]:
        """

        Parameters
        ----------
        header: str

        Returns
        -------
        List[Table]
        """
        return Table().read(self.fits_file, hdu=header)


class ReadoutMATISSE(ReadoutFits):
    """"""
    def __init__(self, *args, **kwargs) -> None:
        """"""
        super().__init__(*args, **kwargs)
        self._oi_flux = None
        self._oi_vis = None
        self._oi_vis2 = None
        self._oi_t3 = None

        self._delay_lines, self._triangles = None, None

    @property
    def oi_flux(self):
        """Fetches the flux"""
        if self._oi_flux is None:
            self._oi_flux = self.get_table_for_fits("oi_flux")
        return self._oi_flux

    @property
    def oi_vis(self):
        """Fetches the visibility"""
        if self._oi_vis is None:
            self._oi_vis = self.get_table_for_fits("oi_vis")
        return self._oi_vis

    @property
    def oi_vis2(self):
        """Fetches the visibility"""
        if self._oi_vis2 is None:
            self._oi_vis2 = self.get_table_for_fits("oi_vis2")
        return self._oi_vis2

    @property
    def oi_t3(self):
        """Fetches the visibility"""
        if self._oi_t3 is None:
            self._oi_t3 = self.get_table_for_fits("oi_t3")
        return self._oi_t3

    @property
    def delay_lines(self):
        """Fetches the delay lines' telescope configuration"""
        if self._delay_lines is None:
            self._delay_lines = ["-".join(list(map(self.telescopes_to_station_names.get,
                                                    station_index)))\
                                  for station_index in self.oi_vis["STA_INDEX"]]
        return self._delay_lines

    @property
    def triangles(self):
        """Fetches the triangles' telescope configurations"""
        if self._triangles is None:
            self._triangles = ["-".join(list(map(self.telescopes_to_station_names.get,
                                                  station_index)))\
                                for station_index in self.oi_t3["STA_INDEX"]]
        return self._triangles

    def get_vis4wl(self, wl_ind: int) -> np.ndarray:
        """Fetches the visdata(amp/phase)/correlated fluxes for one specific wavelength

        Returns
        --------
        visamp4wl: np.ndarray
            The visamp for a specific wavelength
        visamperr4wl: np.ndarray
            The visamperr for a specific wavelength
        visphase4wl: np.ndarray
            The visphase for a specific wavelength
        visphaseerr4wl: np.ndarray
            The visphaseerr for a specific wavelength
        """
        visdata = self.get_vis()
        visamp, visamperr = map(lambda x: x[:6], visdata[:2])
        visphase, visphaseerr = map(lambda x: x[:6], visdata[2:4])
        visamp4wl, visamperr4wl = map(lambda x: np.array([i[wl_ind] for i in x]).flatten(),
                                      [visamp, visamperr])
        visphase4wl, visphaseerr4wl= map(lambda x: np.array([i[wl_ind] for i in x]).flatten(),
                                         [visphase, visphaseerr])

        return visamp4wl, visamperr4wl, visphase4wl, visphaseerr4wl

    def get_vis24wl(self, wl_ind: int) -> np.ndarray:
        """Fetches the vis2data for one specific wavelength

        Returns
        --------
        vis2data4wl: np.ndarray
            The vis2data for a specific wavelength
        vis2err4wl: np.ndarray
            The vis2err for a specific wavelength
        """
        vis2data, vis2err  = map(lambda x: x[:6], self.get_vis2()[:2])
        vis2data4wl, vis2err4wl = map(lambda x: np.array([i[wl_ind] for i in x]).flatten(),
                                      [vis2data, vis2err])

        return vis2data4wl, vis2err4wl


if __name__ == "__main__":
    fits_file = "hd142666_2019-05-14T05_28_03_N_TARGET_FINAL_INT.fits"
    fits_file = DATA_DIR / "jozsef_reductions" / fits_file
    readout = ReadoutMATISSE(fits_file)
    breakpoint()

