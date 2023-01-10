import warnings
from pathlib import Path
from typing import Optional, Dict, List

import numpy as np
import astropy.units as u
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy.coordinates import SkyCoord

# NOTE: Remove units warning. In (.fits)-file 'day' unit contain, which doesn't exist
warnings.simplefilter("ignore", category=u.UnitsWarning)


DATA_DIR = Path(__file__).parent.parent.parent / "data"

# All baseline configurations
SMALL = {1: "A0", 5: "B2", 13: "D0", 10: "C1"}
MED = {28: "K0", 18: "G1", 13: "D0", 24: "J3"}
LARGE = {1: "A0", 18: "G1", 23: "J2", 24: "J3"}
UT = {32: "UT1", 33: "UT2", 34: "UT3", 35: "UT4"}

STATION_INDICES_TO_NAMES = {}
STATION_INDICES_TO_NAMES.update(SMALL)
STATION_INDICES_TO_NAMES.update(MED)
STATION_INDICES_TO_NAMES.update(LARGE)
STATION_INDICES_TO_NAMES.update(UT)


# TODO: Improve docs
# TODO: Add to fluxcalibration that it changes the unit to Jy not ADU -> Jozsef's script
class ReadoutFits:
    """All functionality to work with '.oifits/.fits'-files in a table based manner"""

    def __init__(self, fits_file: Path, flux_file: Optional[Path] = "") -> None:
        self.fits_file = Path(fits_file)
        self.flux_file = Path(flux_file) if flux_file else None

        self._name = None

        self._oi_wl = None
        self._oi_flux, self._oi_t3 = None, None
        self._oi_vis, self._oi_vis2 = None, None
        self._telescope_to_station_names = None

        with fits.open(self.fits_file) as hdul:
            self.primary_header = hdul[0].header

        self.target_name = self.primary_header["OBJECT"]
        self.ra, self.dec = self.primary_header["RA"], self.primary_header["DEC"]
        self.coords = SkyCoord(self.ra*u.deg, self.dec*u.deg, frame="icrs")

    @property
    def name(self):
        """Fetches the target's name via Simbad by its coordinates"""
        # TODO: Make this DISCLAIMER better -> Internet requirement
        if self._name is None:
            header_name = self.primary_header["OBJECT"]
            if header_name in ["SKY", "STD", "STD,RMNREC"]:
                objects = Simbad.query_region(self.coords,
                                              radius=20*u.arcsec)["MAIN_ID"].data.tolist()
                # TODO: Maybe improve the way the common string is found? -> Multiple string names
                self._name = sorted(objects)[0]
            else:
                self._name = header_name.lower()
        return self._name

    @property
    def longest_entry(self):
        """The longest entry of all the rows. Fetched from the 'oi_wl'-tables"""
        return np.max(self.oi_wl["EFF_WAVE"].shape)

    @property
    def oi_wl(self):
        """Gets the wavelength table and reforms it into one entry"""
        if self._oi_wl is None:
            self._oi_wl = Table()
            wl = self.get_table_for_fits("oi_wavelength")["EFF_WAVE"]
            self._oi_wl.add_column(self._oi_wl.Column([wl.data.astype(np.float64)],
                                                      unit=wl.unit), name="EFF_WAVE")
            self._oi_wl["EFF_WAVE"] = self._oi_wl["EFF_WAVE"].to(u.um)
        return self._oi_wl

    @property
    def oi_flux(self):
        """Fetches the flux table if given, and if not makes an empty one"""
        if self._oi_flux is None:
            # NOTE: Not all MATISSE datasets contain 'oi_flux'-data, thus try-except
            try:
                self._oi_flux = self.get_table_for_fits("oi_flux")
            except KeyError:
                # TODO: Implement here flux-file support
                if self.flux_file is not None:
                    ...
                else:
                    self._oi_flux = Table()
                    # TODO: Make this work so the unit is Jy -> Right now it has no effect
                    nan_array = self._oi_flux.Column(np.full(self.longest_entry, np.nan),
                                                     unit=u.Jy)
                    self._oi_flux.add_columns([[nan_array], [nan_array]],
                                              names=["FLUXDATA", "FLUXERR"])
            self._oi_flux.add_column([self.get_table_for_fits("oi_array")["TEL_NAME"].astype(str)],
                                     name="TEL_NAME")
            self._oi_flux.keep_columns(["FLUXDATA", "FLUXERR", "TEL_NAME"])
        return self._oi_flux

    @property
    def oi_vis(self):
        """Fetches the visibility table"""
        if self._oi_vis is None:
            self._oi_vis = self.get_table_for_fits("oi_vis")
            self._oi_vis.add_columns([self.get_delay_lines(self._oi_vis),
                                      self.merge_uv_coords(self._oi_vis),
                                      self.get_baselines(self._oi_vis)],
                                     names=["DELAY_LINE", "UVCOORD", "BASELINE"])
            self._oi_vis.keep_columns(["VISAMP", "VISAMPERR",
                                        "UVCOORD", "DELAY_LINE", "BASELINE"])
        return self._oi_vis

    @property
    def oi_vis2(self):
        """Fetches the squared visibility table"""
        if self._oi_vis2 is None:
            self._oi_vis2 = self.get_table_for_fits("oi_vis2")
            self._oi_vis2.add_columns([self.get_delay_lines(self._oi_vis2),
                                       self.merge_uv_coords(self._oi_vis2),
                                       self.get_baselines(self._oi_vis2)],
                                      names=["DELAY_LINE", "UVCOORD", "BASELINE"])
            self._oi_vis2.keep_columns(["VIS2DATA", "VIS2ERR",
                                        "UVCOORD", "DELAY_LINE", "BASELINE"])
        return self._oi_vis2

    @property
    def oi_t3(self):
        """Fetches the closure phase table"""
        if self._oi_t3 is None:
            self._oi_t3 = self.get_table_for_fits("oi_t3")
            u1, u2 = self._oi_t3["U1COORD"], self._oi_t3["U2COORD"]
            v1, v2 = self._oi_t3["V1COORD"], self._oi_t3["V2COORD"]
            uv_coords = []
            # NOTE: After Jozsef: u3, v3 = -(u1+u2), -(v1+v2) -> Dropping the minus
            # better closure phases in modelling -> Check that!
            for u_coord, v_coord in zip(zip(u1, u2, u1+u2), zip(v1, v2, v1+v2)):
                uv_coords.append(np.array(list(zip(u_coord, v_coord))))
            uv_coords = np.array(uv_coords)
            baselines = [np.sqrt(uv_coord[:, 0]**2+uv_coord[:, 1]**2)\
                    for uv_coord in uv_coords]
            self._oi_t3.add_columns([uv_coords,
                                     self.get_delay_lines(self._oi_t3), baselines],
                                    names=["UVCOORD", "TRIANGLE", "BASELINE"])
            self._oi_t3.keep_columns(["T3PHI", "T3PHIERR",
                                      "UVCOORD", "TRIANGLE", "BASELINE"])
        return self._oi_t3

    @property
    def bcd_configuration(self):
        """Gets the BCD-configuration from the primary header"""
        return ["-".join([self.primary_header["HIERARCH ESO INS BCD1 ID"],
                          self.primary_header["HIERARCH ESO INS BCD2 ID"]])]

    @property
    def tpl_start(self):
        """Gets the template start datetime from the primary header"""
        return self.primary_header["HIERARCH ESO TPL START"]

    def is_calibrator(self):
        """Checks if the target was observed as calibrator with info from the primary
        header"""
        if self.primary_header["HIERARCH ESO DPR CATG"].lower() == "calib":
            return True
        return False

    def get_table_for_fits(self, header: str) -> List[Table]:
        """Fetches the data from a (.fits)-file as a table

        Parameters
        ----------
        header: str
            The header which data is to be stored in a table

        Returns
        -------
        Table
        """
        return Table().read(self.fits_file, hdu=header)

    def merge_uv_coords(self, table: Table):
        """Merges the u- and v-coordinates into a set of (u, v)-coordinates"""
        u_coords, v_coords = table["UCOORD"], table["VCOORD"]
        return np.array(list(zip(u_coords, v_coords)))

    def get_baselines(self, table: Table):
        """Calculates the baselines from the (u, v)-coordinates"""
        u_coords, v_coords = table["UCOORD"], table["VCOORD"]
        return np.sqrt(u_coords**2+v_coords**2)

    def get_delay_lines(self, table: Table):
        """Fetches the delay lines' telescope configuration from the visibility table"""
        return ["-".join(list(map(STATION_INDICES_TO_NAMES.get, station_index)))
                for station_index in table["STA_INDEX"]]


if __name__ == "__main__":
    fits_files = ["HD_163296_2019-03-23T08_41_19_N_TARGET_FINALCAL_INT.fits",
                  "HD_163296_2019-03-23T08_41_19_L_TARGET_FINALCAL_INT.fits",
                  "HD_163296_2019-05-06T08_19_51_L_TARGET_FINALCAL_INT.fits"]
    fits_files = [DATA_DIR / "tests" / fits_file for fits_file in fits_files]
    readout = ReadoutFits(fits_files[0])
    breakpoint()
