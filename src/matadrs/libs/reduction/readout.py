from typing import Dict, List

import numpy as np
import astropy.units as u
from astropy.io import fits
from pathlib import Path
from astropy.table import Table, vstack


DATA_DIR = Path(__file__).parent.parent.parent.parent.parent / "data"


class ReadoutFits:
    """All functionality to work with '.oifits/.fits'-files"""
    def __init__(self, fits_files: List[Path]) -> None:
        self.hdulists = None
        self._telescope_to_station_names = None

        if not isinstance(fits_files, list):
            fits_files = [fits_files]

        if not isinstance(fits_files[0], Path):
            self.fits_files = [Path(fits_file) for fits_file in fits_files]
        else:
            self.fits_files = fits_files

    def __enter__(self) -> None:
        """"""
        self.hdulists = [fits.open(fits_file) for fits_file in self.fits_files]
        return self

    def __exit__(self, *exc) -> None:
        """"""
        for hdul in self.hdulist:
            hdul.close()

    def __repr__(self):
        """The class's representation"""
        # TODO: Make this better
        return f"The instance is maintaining the following files:\n"\
               f"{', '.join([fits_file.name for fits_file in self.fits_files])}"

    def __str__(self):
        """The class's print"""
        return self.__repr__()

    @property
    def telescopes_to_station_names(self) -> List[Dict]:
        """"""
        if self._telescope_to_station_names is None:
            tel_names = self.get_tables_for_fits("oi_array")
            self._telescope_to_station_names = ...
        return self._telescope_to_station_names

    def get_primary_header(self):
        """"""
        return [hdul[0].header for hdul in self.hdulists]

    def get_target(self):
        """"""
        return [primary["OBJECT"] for primary in self.get_primary_header()]

    def get_ra_and_dec(self):
        """"""
        return [(primary["RA"], primary["DEC"]) for primary in self.get_primary_header()]

    def get_bcd_configuration(self):
        """"""
        return ["-".join([primary["HIERARCH ESO INS BCD1 ID"],
                          primary["HIERARCH ESO INS BCD2 ID"]])\
                for primary in self.hdulists]

    # def get_split_uvcoords(self) -> np.ndarray:
        # """Splits a 2D-np.array into its 1D-components and returns the u- and
        # v-coords seperatly"""
        # uvcoords = self.get_uvcoords()
        # return np.array([item[0] for item in uvcoords]),\
                # np.array([item[1] for item in uvcoords])

    def get_tables_for_fits(self, header: str) -> List[Table]:
        """

        Parameters
        ----------
        header: str

        Returns
        -------
        List[Table]
        """
        return [Table().read(fits_file, hdu=header) for fits_file in self.fits_files]


class ReadoutMATISSE(ReadoutFits):
    """"""
    def __init__(self, *args, **kwargs) -> None:
        """"""
        super().__init__(*args, **kwargs)
        self._oi_vis = None

    @property
    def oi_flux(self):
        """Fetches the flux"""
        if self._oi_flux is None:
            self._oi_flux = self.get_tables_for_fits("oi_flux")
        return self._oi_flux

    @property
    def oi_vis(self):
        """Fetches the visibility"""
        if self._oi_vis is None:
            self._oi_vis = self.get_tables_for_fits("oi_vis")
        return self._oi_vis

    @property
    def oi_vis2(self):
        """Fetches the visibility"""
        if self._oi_vis2 is None:
            self._oi_vis2 = self.get_tables_for_fits("oi_vis2")
        return self._oi_vis2

    @property
    def oi_t3(self):
        """Fetches the visibility"""
        if self._oi_t3 is None:
            self._oi_t3 = self.get_tables_for_fits("oi_t3")
        return self._oi_t3

    @property
    def wl(self):
        """Fetches the wavelengths"""
        if self._wl is None:
            self._wl = self.get_tables_for_fits("oi_wl")["EFF_WAVE"].to(u.um)
        return self._wl

    # self.tel_vis = [["-".join(list(map(TEL_STATIONS.get, station_index["STA_INDEX"])))\
                    # for station_index in vis] for vis in self.oi_vis]

    # self.tel_t3 = [["-".join(list(map(TEL_STATIONS.get, station_index["STA_INDEX"])))\
                    # for station_index in t3] for t3 in self.oi_t3]

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
    fits_files = ["hd142666_2019-05-14T05_28_03_N_TARGET_FINAL_INT.fits",
                  "hd163296_2019-03-23T08_41_19_L_TARGET_FINAL_INT.fits"]
    fits_files = [DATA_DIR / "jozsef_reductions" / fits_file for fits_file in fits_files]
    with ReadoutMATISSE(fits_files) as readout:
        breakpoint()

