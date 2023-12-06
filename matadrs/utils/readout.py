import pkg_resources
import re
import warnings
from logging import warn
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import astropy.units as u
from astropy.io import fits
from astroquery.simbad import Simbad
from astropy.table import Table
from astropy.coordinates import SkyCoord

# NOTE: Remove units warning. In (.fits)-file 'day' unit contain, which doesn't exist
warnings.simplefilter("ignore", category=u.UnitsWarning)

__all__ = ["ReadoutFits"]

DATA_DIR = Path(pkg_resources.resource_filename("matadrs", "data"))
GRAVITY_TO_INDEX = {"sc": 10, "ft": 20}


class ReadoutFits:
    """Reads out the Cards as Tables as well as the primary header of a
    (.fits)-file and makes certain keys from the primary header available
    as properties of the class.

    Parameters
    ----------
    fits_file : pathlib.Path
        The (.fits)-file from which data is sourced.

    Attributes
    ----------
    fits_file : pathlib.Path
        The (.fits)-file that has been read in.
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

    def __init__(self, fits_file: Path) -> None:
        """The class's constructor"""
        self.fits_file = Path(fits_file)

        self._sta_to_tel = None
        self._name, self._coords = None, None
        self._simbad_query = None
        self.gravity_method = "ft"

        headers = ["oi_array", "oi_wavelength", "oi_flux",
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
                if "OBJECT" in self.primary_header:
                    header_name = self.primary_header["OBJECT"]
                else:
                    header_name = self.primary_header["HIERARCH ESO OBS TARG NAME"]

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
    def gravity_index(self) -> int:
        """Returns the indices for either the fringe tracker or science
        observations."""
        if self.primary_header["instrume"].lower() == "gravity":
            return GRAVITY_TO_INDEX[self.gravity_method]
        return None

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
    def pipeline_version(self) -> str:
        """Fetches the pipeline version from the primary header."""
        return self.primary_header["HIERARCH ESO PRO REC1 PIPE ID"]

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
    def resolution(self) -> str:
        """Fetches the object's N-band resolutions from the primary header."""
        return self.primary_header["HIERARCH ESO INS DIN NAME"].lower()

    @property
    def longest_entry(self) -> int:
        """The longest entry of all the rows fetched from the 'oi_wl' Table."""
        return np.max(self.oi_wavelength["EFF_WAVE"].shape)

    @property
    def sta_to_tel(self) -> Dict[int, str]:
        """Gets the telescope's station index to telescope name mapping."""
        if self._sta_to_tel is None:
            self._sta_to_tel = dict(zip(self.oi_array["STA_INDEX"],
                                        self.oi_array["STA_NAME"]))
        return self._sta_to_tel

    @property
    def oi_wavelength(self) -> Table:
        """Gets the wavelength table and reforms it into one entry."""
        if self._oi_wavelength is None:
            self._oi_wavelength = Table()
            wl = self.get_table_for_fits("oi_wavelength")["EFF_WAVE"]
            self._oi_wavelength.add_column(self._oi_wavelength.Column(
                [wl.data.astype(np.float64)], unit=wl.unit), name="EFF_WAVE")
            if self.oi_wavelength["EFF_WAVE"].unit is not u.m:
                self.oi_wavelength["EFF_WAVE"] = self.oi_wavelength["EFF_WAVE"].value*u.m
            self._oi_wavelength["EFF_WAVE"] = self._oi_wavelength["EFF_WAVE"].to(u.um)
        return self._oi_wavelength

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
            # NOTE: Not all MATISSE datasets contain 'oi_flux'-data.
            # Thus try-except
            try:
                self._oi_flux = self.get_table_for_fits("oi_flux")
            except KeyError:
                self._oi_flux = Table()
                nan_array = self._oi_flux.Column(
                        np.full(self.longest_entry, np.nan))
                nan_array.unit = u.Jy
                self._oi_flux.add_columns(
                        [[nan_array], [nan_array], [np.nan]],
                        names=["FLUXDATA", "FLUXERR", "STA_INDEX"])
            if "FLUXDATA" in self._oi_flux.columns:
                self._oi_flux.keep_columns(
                        ["FLUXDATA", "FLUXERR", "STA_INDEX"])
            else:
                self._oi_flux.keep_columns(
                        ["FLUX", "FLUXERR", "STA_INDEX"])
        return self._oi_flux

    @property
    def oi_vis(self) -> Table:
        """Fetches the visibility table."""
        if self._oi_vis is None:
            try:
                self._oi_vis = self.get_table_for_fits("oi_vis")
                self._oi_vis.add_columns(
                        [self.get_delay_lines(self._oi_vis),
                         self.merge_uv_coords(self._oi_vis),
                         self.get_baselines(self._oi_vis)],
                        names=["DELAY_LINE", "UVCOORD", "BASELINE"])
            except KeyError:
                self._oi_vis = Table()
                nan_column = self._oi_vis.Column(
                        np.full(self.longest_entry, np.nan))
                nan_array = [nan_column for _ in range(6)]
                nan_six = [(np.nan, np.nan) for _ in range(6)]
                nan_str = ["" for _ in range(6)]
                nan_base = [np.nan for _ in range(6)]
                self._oi_vis.add_columns(
                        [nan_array, nan_array, nan_six,
                         nan_array, nan_array, nan_str,
                         nan_base, nan_base, nan_array, nan_six],
                        names=["VISAMP", "VISAMPERR", "UVCOORD",
                               "VISPHI", "VISPHIERR",
                               "DELAY_LINE", "BASELINE",
                               "MJD", "FLAG", "STA_INDEX"])
            self._oi_vis.keep_columns(
                    ["VISAMP", "VISAMPERR", "UVCOORD",
                     "VISPHI", "VISPHIERR",
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
        return self.observation_type == "calib"

    def is_pip_version_greater_equal(self, version: str) -> bool:
        """Checks if the pipeline's version is greater than the
        reference version.

        Parameters
        ----------
        version : str
            The version to compare with the reference version.
            To be passed in the format 'x.y.z'

        Returns
        -------
        greater : bool
        """
        numbers_to_check = [int(num) for num
                            in re.findall(r'\d+', version)]
        numbers_reference = [int(num) for num
                             in re.findall(r'\d+', self.pipeline_version)]

        if len(numbers_to_check) >= 3 and len(numbers_reference) >= 3:
            return all(numbers_reference[index] >= numbers_to_check[index]
                       for index in range(3))
        else:
            raise ValueError("Invalid version format."
                             " Please use x.y.z format.")

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
        unit = getattr(self, header)[sub_header.upper()].unit
        return str(unit) if unit is not None else "a.u."

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
        if header == "oi_array":
            return Table().read(self.fits_file, hdu=header)
        with fits.open(self.fits_file, "readonly") as hdul:
            return Table().read(hdul[header, self.gravity_index])

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

    def get_telescopes(self) -> List[str]:
        """Fetches the station indices from the 'oi_flux' table
        and returns the telescope's names.

        Returns
        -------
        telescopes : list of  str
        """
        return [self.sta_to_tel[tel] for tel in self.oi_flux["STA_INDEX"]]

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
