from logging import warn
import pkg_resources
import warnings
from pathlib import Path
from typing import Dict, Tuple, List, Optional

import numpy as np
import astropy.units as u
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.interpolate import CubicSpline

# NOTE: Remove units warning. In (.fits)-file 'day' unit contain, which doesn't exist
warnings.simplefilter("ignore", category=u.UnitsWarning)

__all__ = ["ReadoutFits"]

DATA_DIR = Path(pkg_resources.resource_filename("matadrs", "data"))


# TODO: Add to fluxcalibration that it changes the unit to Jy not ADU -> Maybe in Jozsef's
# script?
class ReadoutFits:
    """Reads out the Cards as Tables as well as the primary header of a
    (.fits)-file and makes certain keys from the primary header available
    as properties of the class.

    Parameters
    ----------
    fits_file : pathlib.Path
        The (.fits)-file from which data is sourced.
    flux_file : pathlib.Path, optional
        Additional flux file, that replaces the total flux data for the
        (.fits)-file.

    Attributes
    ----------
    fits_file : pathlib.Path
        The (.fits)-file that has been read in.
    flux_file : pathlib.Path, optional
        If provided, substitutes the columns of the 'oi_flux' Table with the
        values from the flux file.
    primary_header : astropy.io.fits.PrimaryHDU
        The primary header of the (.fits)-file.
    name
    ra
    dec
    coords
    observation_type
    array_configuration
    bcd_configuration
    tpl_start
    resolution
    detector
    seeing
    tau0
    longest_entry
    oi_wl
    oi_wl_hdr
    oi_array
    oi_array_hdr
    oi_flux
    oi_flux_hdr
    oi_t3
    oi_t3_hdr
    oi_vis
    oi_vis_hdr
    oi_vis2
    oi_vis2_hdr

    Methods
    -------
    is_calibrator()
        Fetches the object's observation mode and returns true if it has been
        observed in 'CALIB' mode.
    get_flux_data_from_flux_file():
        Reads the flux data from the flux file and then interpolates it to the
        wavelength solution used by MATISSE.
    get_table_for_fits(header: str)
        Fetches a Card by its header and then reads its information into a
        Table.
    merge_uv_coords(table)
        Fetches the u- and v-coordinates from a Table and the merges them into
        a set of (u, v)-coordinates.
    get_baselines(table)
        Fetches the u- and v-coordinates from a Table and calculates their
        baselines.
    get_delay_lines(table)
        Fetches the station indices from a Table and returns the telescope's
        delay line configuration.
    """

    def __init__(self, fits_file: Path, flux_file: Optional[Path] = "") -> None:
        """The class's constructor"""
        # TODO: Maybe rename fits file to file? or Path?
        self.fits_file = Path(fits_file)
        self.flux_file = Path(flux_file) if flux_file else None

        self._sta_to_tel = None
        self._name, self._coords = None, None
        self._simbad_query = None

        headers = ["oi_array", "oi_wl", "oi_flux",
                   "oi_t3", "oi_vis", "oi_vis2"]
        for header in headers:
            setattr(self, f"_{header}", None)
            setattr(self, f"{header}_hdr", self.get_header(header))

        with fits.open(self.fits_file) as hdul:
            self.primary_header = hdul[0].header

    @property
    def name(self) -> str:
        """Fetches the object's name from the primary header and if not found or not named
        in it tries it via Simbad by its coordinates.

        Notes
        -----
        Fetching the name via simbad by its coordinates REQUIRES online access.
        """
        if self._name is None:
            try:
                header_name = self.primary_header["OBJECT"]
            except KeyError:
                header_name = None
            if (header_name in ["SKY", "STD", "STD,RMNREC"]) or (header_name is None):
                objects = Simbad.query_region(self.coords,
                                              radius=20*u.arcsec)["MAIN_ID"].data.tolist()
                self._name = sorted(objects)[0]
            else:
                self._name = header_name.lower()
        return self._name

    @property
    def simbad_query(self):
        """The simbad_query property."""
        if self._simbad_query is None:
            simbad = Simbad.query_object(self.name)
            ra, dec = simbad["RA"], simbad["DEC"]
            self._simbad_query = SkyCoord(ra, dec, unit="deg")
        return self._simbad_query

    @property
    def ra(self) -> str:
        """Fetches the right ascension from the primary header."""
        if "RA" in self.primary_header:
            return self.primary_header["RA"]
        else:
            return self.simbad_query.ra

    @property
    def dec(self) -> str:
        """Fetches the declination from the primary header."""
        if "DEC" in self.primary_header:
            return self.primary_header["DEC"]
        else:
            return self.simbad_query.dec

    @property
    def mjd(self) -> str:
        """Fetches the observation's modified julian date from the primary header."""
        if "MJD-OBS" in self.primary_header:
            return self.primary_header["MJD-OBS"]
        return None

    @property
    def coords(self) -> SkyCoord:
        """Fetches both right ascension and declination from the primary header and wraps
        it via astropy's Skycoord class."""
        if self._coords is None:
            self._coords = SkyCoord(self.ra*u.deg, self.dec*u.deg, frame="icrs")
        return self._coords

    @property
    def observation_type(self) -> str:
        """Fetches the type of the observation, i.e., if the object is a science target or
        calibrator."""
        if "SCI" in self.primary_header["HIERARCH ESO OBS NAME"]:
            return "science"
        else:
            return self.primary_header["HIERARCH ESO DPR CATG"].lower()

    @property
    def array_configuration(self) -> str:
        """Fetches the array's configuration from the primary header."""
        if "HIERARCH ESO ISS BASELINE" in self.primary_header:
            array = self.primary_header["HIERARCH ESO ISS BASELINE"]
        else:
            array = self.primary_header["HIERARCH ESO OBS BASELINE"]
        return "uts" if "UT" in array else "ats"

    @property
    def bcd_configuration(self) -> str:
        """Fetches the BCD-configuration from the primary header."""
        return "-".join([self.primary_header["HIERARCH ESO INS BCD1 ID"],
                         self.primary_header["HIERARCH ESO INS BCD2 ID"]]).lower()

    @property
    def tpl_start(self) -> str:
        """Fetches the observation's start datetime from the primary header."""
        return self.primary_header["HIERARCH ESO TPL START"]

    @property
    def detector(self) -> str:
        """Fetches the detector used for the observation from the primary header."""
        return self.primary_header["HIERARCH ESO DET NAME"]

    @property
    def seeing(self):
        """Fetches the seeing from the primary header."""
        return self.primary_header["HIERARCH ESO ISS AMBI FWHM START"]

    @property
    def tau0(self):
        """Fetches the tau0 from the primary header."""
        return self.primary_header["HIERARCH ESO ISS AMBI TAU0 START"]

    @property
    def resolution(self) -> Tuple[str, str]:
        """Fetches the object's N-band resolutions from the primary header."""
        return self.primary_header["HIERARCH ESO INS DIN NAME"].lower()

    @property
    def longest_entry(self) -> int:
        """The longest entry of all the rows fetched from the 'oi_wl' Table."""
        return np.max(self.oi_wl["EFF_WAVE"].shape)

    @property
    def sta_to_tel(self) -> Dict[int, str]:
        """Gets the telescope's station index to telescope name mapping."""
        if self._sta_to_tel is None:
            self._sta_to_tel = dict(zip(self.oi_array["STA_INDEX"],
                                        self.oi_array["STA_NAME"]))
        return self._sta_to_tel

    @property
    def oi_wl(self) -> Table:
        """Gets the wavelength table and reforms it into one entry."""
        if self._oi_wl is None:
            self._oi_wl = Table()
            wl = self.get_table_for_fits("oi_wavelength")["EFF_WAVE"]
            self._oi_wl.add_column(self._oi_wl.Column([wl.data.astype(np.float64)],
                                                      unit=wl.unit), name="EFF_WAVE")
            if self.oi_wl["EFF_WAVE"].unit is not u.m:
                self.oi_wl["EFF_WAVE"] = self.oi_wl["EFF_WAVE"].value*u.m
            self._oi_wl["EFF_WAVE"] = self._oi_wl["EFF_WAVE"].to(u.um)
        return self._oi_wl

    @property
    def oi_array(self) -> Table:
        """Fetches the array's information."""
        if self._oi_array is None:
            self._oi_array = self.get_table_for_fits("oi_array")
        return self._oi_array

    @property
    def oi_flux(self) -> Table:
        """Fetches the flux table if given, and if not makes an empty one."""
        if self._oi_flux is None:
            # NOTE: Not all MATISSE datasets contain 'oi_flux'-data, thus try-except
            try:
                self._oi_flux = self.get_table_for_fits("oi_flux")
            except KeyError:
                self._oi_flux = Table()
                if self.flux_file not in ["", None]:
                    # TODO: Make this work so the unit is Jy -> Right now it has no effect
                    nan_array = self._oi_flux.Column(np.full(self.longest_entry, np.nan),
                                                     unit=u.Jy)
                    self._oi_flux.add_columns([[nan_array], [nan_array]],
                                              names=["FLUXDATA", "FLUXERR"])
                else:
                    flux, flux_err = self.get_flux_data_from_flux_file()
                    self._oi_flux.add_columns([self._oi_flux.Column([flux], unit=u.Jy),
                                              self._oi_flux.Column([flux_err], unit=u.Jy)],
                                              names=["FLUXDATA", "FLUXERR"])
            self._oi_flux.keep_columns(["FLUXDATA", "FLUXERR"])
        return self._oi_flux

    @property
    def oi_vis(self) -> Table:
        """Fetches the visibility table."""
        if self._oi_vis is None:
            self._oi_vis = self.get_table_for_fits("oi_vis")
            self._oi_vis.add_columns([self.get_delay_lines(self._oi_vis),
                                      self.merge_uv_coords(self._oi_vis),
                                      self.get_baselines(self._oi_vis)],
                                     names=["DELAY_LINE", "UVCOORD", "BASELINE"])
            self._oi_vis.keep_columns(["VISAMP", "VISAMPERR", "UVCOORD",
                                       "DELAY_LINE", "BASELINE",
                                       "MJD", "FLAG", "STA_INDEX"])
        return self._oi_vis

    @property
    def oi_vis2(self) -> Table:
        """Fetches the squared visibility table."""
        if self._oi_vis2 is None:
            self._oi_vis2 = self.get_table_for_fits("oi_vis2")
            self._oi_vis2.add_columns([self.get_delay_lines(self._oi_vis2),
                                       self.merge_uv_coords(self._oi_vis2),
                                       self.get_baselines(self._oi_vis2)],
                                      names=["DELAY_LINE", "UVCOORD", "BASELINE"])
            self._oi_vis2.keep_columns(["VIS2DATA", "VIS2ERR",
                                        "UVCOORD", "DELAY_LINE",
                                        "BASELINE", "MJD", "FLAG", "STA_INDEX"])
        return self._oi_vis2

    @property
    def oi_t3(self) -> Table:
        """Fetches the closure phase table."""
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

    def is_calibrator(self) -> bool:
        """Fetches the object's observation mode and returns true if it has been observed
        in 'CALIB' mode.

        Returns
        -------
        observed_as_calibrator : bool
        """
        if self.observation_type == "calib":
            return True
        return False

    # TODO: Get a better error representation for the flux
    def get_flux_data_from_flux_file(self) -> Tuple[u.Quantity[u.Jy],
                                                    u.Quantity[u.Jy]]:
        """Reads the flux data from the flux file and then interpolates it to the
        wavelength solution used by MATISSE.

        Returns
        -------
        flux : astropy.units.Jy
            The, to MATISSE's wavelength solution interpolated, flux fetched from the
            flux file.
        flux_error : astropy.units.Jy
            The, to MATISSE's wavelength solution interpolated, flux error fetched from
            the flux file.
        """
        flux_data = Table.read(self.flux_file, names=["wl", "flux"], format="ascii")
        cubic_spline = CubicSpline(flux_data["wl"], flux_data["flux"])
        interpolated_flux = (cubic_spline(self.oi_wl["EFF_WAVE"].data)).ravel()
        return interpolated_flux, interpolated_flux*0.1

    def get_header(self, header: str) -> Optional[fits.Header]:
        """Fetches a Card's header by its header name.

        Parameters
        ----------
        header : str
            The header which data is to be stored in a table.

        Returns
        -------
        header : fits.Header, optional
        """
        with fits.open(self.fits_file, "readonly") as hdul:
            if header not in hdul:
                warn(f"Header {header} not found!")
                return
            return hdul[header].header

    def get_unit(self, header: str, sub_header: str) -> str:
        """Fetches the unit of a header by the sub header's name."""
        if header == "oi_vis"\
                and header["AMPTYP"] == "correlated flux":
            return "Jy"
        header = self.get_header(header)
        list_index = list(header.values()).index(sub_header.upper())
        key = f"TUNIT{list(header.keys())[list_index][-1]}"
        if key in header:
            return header[key]
        return "a.u."

    def get_table_for_fits(self, header: str) -> Table:
        """Fetches a Card by its header and then reads its information into a Table.

        Parameters
        ----------
        header : str
            The header which data is to be stored in a table.

        Returns
        -------
        Table
        """
        return Table().read(self.fits_file, hdu=header)

    def merge_uv_coords(self, table: Table) -> np.ndarray:
        """Fetches the u- and v-coordinates from a Table and the merges them into a set
        of (u, v)-coordinates.

        Parameters
        ----------
        table : Table
            The Table to be read in.

        Returns
        -------
        merged_uv_coords : numpy.ndarray
        """
        return np.array(list(zip(table["UCOORD"], table["VCOORD"])))

    def get_baselines(self, table: Table) -> np.ndarray:
        """Fetches the u- and v-coordinates from a Table and calculates their baselines.

        Parameters
        ----------
        table : Table
            The Table to be read in.

        Returns
        -------
        baselines  : numpy.ndarray
        """
        return np.sqrt(table["UCOORD"]**2+table["VCOORD"]**2)

    def get_delay_lines(self, table: Table) -> List[str]:
        """Fetches the station indices from a Table and returns the telescope's delay
        line configuration.

        Parameters
        ----------
        table : Table
            The Table to be read in.

        Returns
        -------
        delay_lines : list of  str
        """
        return ["-".join(list(map(self.sta_to_tel.get, station_index)))
                if all([index in self.sta_to_tel for index in station_index]) else ""
                for station_index in table["STA_INDEX"]]


if __name__ == "__main__":
    fits_files = ["HD_163296_2019-03-23T08_41_19_N_TARGET_FINALCAL_INT.fits",
                  "HD_163296_2019-03-23T08_41_19_L_TARGET_FINALCAL_INT.fits",
                  "HD_163296_2019-05-06T08_19_51_L_TARGET_FINALCAL_INT.fits"]
    fits_files = [DATA_DIR / "tests" / fits_file for fits_file in fits_files]
    flux_files = ["HD_163296_sws.txt", "HD_163296_timmi2.txt"]
    flux_files = [DATA_DIR / "tests" / flux_file for flux_file in flux_files]
    readout = ReadoutFits(fits_files[0], flux_files[0])
