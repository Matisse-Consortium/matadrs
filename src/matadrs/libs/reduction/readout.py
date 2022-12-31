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

        self._oi_wl = None
        self._oi_flux = None
        self._oi_vis = None
        self._oi_vis2 = None
        self._oi_t3 = None
        self._telescope_to_station_names = None

        with fits.open(self.fits_file) as hdul:
            self.primary_header = hdul[0].header

        self.target_name = self.primary_header["OBJECT"]
        self.ra, self.dec = self.primary_header["RA"], self.primary_header["DEC"]

    @property
    def oi_wl(self):
        """"""
        if self._oi_wl is None:
            self._oi_wl = self.get_table_for_fits("oi_wavelength")
            self._oi_wl.remove_column("EFF_BAND")
        return self._oi_wl

    @property
    def oi_flux(self):
        """Fetches the flux"""
        if self._oi_flux is None:
            self._oi_flux = self.get_table_for_fits("oi_flux")
            self._oi_flux.add_columns([self.get_table_for_fits("oi_array")["TEL_NAME"]],
                                     names="TEL_NAME")
            self._oi_flux.keep_columns(["FLUXDATA", "FLUXERR", "TEL_NAME"])
        return self._oi_flux

    @property
    def oi_vis(self):
        """Fetches the visibility"""
        if self._oi_vis is None:
            self._oi_vis = self.get_table_for_fits("oi_vis")
            self._oi_vis.add_columns([self.get_delay_lines(self._oi_vis)],
                                     names=["DELAY_LINE"])
            self._oi_vis.keep_columns(["VISAMP", "VISAMPERR",
                                       "UCOORD", "VCOORD", "DELAY_LINE"])
        return self._oi_vis

    @property
    def oi_vis2(self):
        """Fetches the visibility"""
        if self._oi_vis2 is None:
            self._oi_vis2 = self.get_table_for_fits("oi_vis2")
            self._oi_vis2.add_columns([self.get_delay_lines(self._oi_vis2)],
                                      names=["DELAY_LINE"])
            self._oi_vis2.keep_columns(["VIS2DATA", "VIS2ERR",
                                        "UCOORD", "VCOORD", "DELAY_LINE"])
        return self._oi_vis2

    @property
    def oi_t3(self):
        """Fetches the closure phase"""
        if self._oi_t3 is None:
            self._oi_t3 = self.get_table_for_fits("oi_t3")
            # NOTE: After Jozsef this does not make good closure phases
            # u3, v3 = -(u1+u2), -(v1+v2)
            u1, u2 = self._oi_t3["U1COORD"].data, self._oi_t3["U2COORD"].data
            v1, v2 = self._oi_t3["V1COORD"].data, self._oi_t3["V2COORD"].data
            self._oi_t3.add_columns([list(zip(u1, u2, u1+u2)),
                                    list(zip(v1, v2, v1+v2)),
                                    self.get_triangles(self._oi_t3)],
                                   names=["UCOORD", "VCOORD", "TRIANGLE"])
            self._oi_t3.keep_columns(["T3PHI", "T3PHIERR",
                                      "UCOORD", "VCOORD", "TRIANGLE"])
        return self._oi_t3

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

    def get_delay_lines(self, oi_vis: Table):
        """Fetches the delay lines' telescope configuration"""
        return ["-".join(list(map(self.telescopes_to_station_names.get,
                                  station_index)))\
                for station_index in oi_vis["STA_INDEX"]]

    def get_triangles(self, oi_t3: Table):
        """Fetches the triangles' telescope configurations"""
        return ["-".join(list(map(self.telescopes_to_station_names.get,
                                  station_index)))\
                for station_index in oi_t3["STA_INDEX"]]


class DataPrep:
    def __init__(self, fits_files: List[Path]) -> None:
        self.fits_files = fits_files


if __name__ == "__main__":
    fits_file = "hd142666_2019-05-14T05_28_03_N_TARGET_FINAL_INT.fits"
    fits_file = DATA_DIR / "jozsef_reductions" / fits_file
    readout = ReadoutMATISSE(fits_file)

