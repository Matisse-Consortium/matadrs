from pathlib import Path
from typing import Callable, Tuple, List, Union, Optional

import astropy.units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation
from astropy.table import Column
from matplotlib.axes import Axes
from pandas import DataFrame

from ..utils.readout import ReadoutFits


# TODO: Make this go away from class and just individual plotting functions? -> Maybe
# easier for paper relevant plots, but first make unified class that makes nice data
# reduction plots
# TODO: Implement text plotter with the information on the observation
class Plotter:
    """Class that plots models as well as reduced data

    Parameters
    ----------
    fits_files : pathlib.Path or list of pathlib.Path
    flux_files : pathlib.Path or list of pathlib.Path, optional
    plot_name : str, optional
    save_path : pathlib.Path, optional

    Attributes
    ----------
    num_components : int

    Methods
    -------
    band_mask(wl)
        The masking for the bands to make the plots visually readable.
    mask_dataframe(df, column_mask)
    calculate_uv_points(baselines, hour_angle, latitude)
    make_uv_tracks(ax, uv_coord, baselines, sta_label, flag,
                   symbol, color, sel_wl, airmass_lim)
    plot_uv(ax, symbol, color, sel_wl, airmass_lim)
    set_dataframe(self, labels, column)
        Prepares each row in a column as a pandas DataFrame.
    make_component(data_name, legend_format, unwrap, period)
        Generates a pandas DataFrame that has all the plots' information.
    add_flux()
        Adds the total flux(es) as a subplot.
    add_vis(corr_flux=False, legend_format="long")
        Adds the visibilities/correlated fluxes as a subplot.
    add_vis2(legend_format="long")
        Adds the squared visibilities as a subplot.
    add_cphases(unwrap=False, period=360)
        Adds the closure phases as a subplot.
    add_diff_phases(unwrap=False, period=360)
        Adds the differential phases as a subplot.
    add_uv()
        Adds the (u, v)-coordinates as a subplot.
    add_text()
        To be implemented.
    add_mosaic(unwrap=False, legend_format="long")
        Combines multiple subplots to produce a mosaic plot.
    get_plot_linestyle(already_chosen_linestyles)
        Gets a linestyle, which is different from the already chosen one.
    plot(self, save=False)
        Combines the individual components into one plot.
    """

    def __init__(self, fits_files: Path | List[Path],
                 flux_files: Optional[Path | List[Path]] = None,
                 plot_name: Optional[str] = None,
                 save_path: Optional[Path] = None) -> None:
        """The class's constructor"""
        self.fits_files = [fits_files]\
            if not isinstance(fits_files, List) else fits_files

        if flux_files is None:
            self.flux_files = [None]*len(self.fits_files)
        elif isinstance(flux_files, List) and len(flux_files) != len(self.fits_files):
            raise IOError("Flux files must either be None or be the same number"
                          "the (.fits)-files provided")
        else:
            self.flux_files = flux_files

        if save_path is None:
            self.save_path = Path("").cwd()
        else:
            self.save_path = Path(save_path)

        if plot_name is not None:
            self.plot_name = plot_name
        elif isinstance(fits_files, List) and len(fits_files) > 1:
            self.plot_name = "combined_fits.pdf"
        else:
            self.plot_name = f"{Path(fits_files).stem}.pdf"

        self.readouts = [ReadoutFits(fits_file, flux_file)
                         for fits_file, flux_file
                         in zip(self.fits_files, self.flux_files)]

        # TODO: Maybe add more linestyles and colors or make this implementation better
        self.linestyles = ["-", "--", "-.", ":"]
        self.colors = ["b", "g", "r", "c", "m", "y"]
        self.components = {}

    def __repr__(self):
        """The class's console representation"""
        return f"Plotting the following (.fits)-files:\n{'':-^50}\n" + \
            "\n".join([readout.fits_file.stem for readout in self.readouts])

    def __str__(self):
        """The class's string representation"""
        return self.__repr__()

    @property
    def num_components(self):
        """The number of componets contained"""
        return len(self.components)

    # TODO: Rename function at a future point
    def calculate_uv_points(self, baselines: List[float],
                            hour_angle: np.ndarray[float],
                            latitude: u.rad,
                            declination: u.rad) -> Tuple[np.ndarray]:
        """Calculates the earth rotation (synthesis) for the uv-point
        corresponding to the baselines for the input hour angle(s)

        Parameters
        -----------
        baselines : list of float
            The baselines in the following order: Baselines east, -north, -longest.
        hour_angle : numpy.ndarray of float
        latitude : astropy.units.rad
            The latitude of the site
        declination : astropy.units.rad

        Returns
        -------
        u_coords : numpy.ndarray
        v_coords : numpy.ndarray
        """
        baseline_east, baseline_north, baseline_longest = baselines

        u_coords = baseline_east * np.cos(hour_angle) - baseline_north * np.sin(latitude)\
            * np.sin(hour_angle) + baseline_longest * np.cos(latitude) * np.sin(hour_angle)
        v_coords = baseline_east * np.sin(declination) * np.sin(hour_angle)\
            + baseline_north * (np.sin(latitude) * np.sin(declination) * np.cos(hour_angle)
                                + np.cos(latitude) * np.cos(declination)) - baseline_longest * \
            (np.cos(latitude) * np.sin(declination) * np.cos(hour_angle)
             - np.sin(latitude) * np.cos(declination))
        return u_coords, v_coords

    # TODO: Add proper docs
    def make_uv_tracks(self, ax, uv_coord: np.ndarray[float],
                       baselines: List[np.ndarray],
                       sta_label: List[np.ndarray],
                       declination: float, flag: bool,
                       symbol: str, color: str, sel_wl: float,
                       airmass_lim: float) -> None:
        """This function was written by Jozsef Varga (from menEWS: menEWS_plot.py).

        From coordinate + ha (range), calculate uv tracks

        Parameters
        ----------
        uv_coords : numpy.ndarray of float
        baselines : list of numpy.ndarray
            The baselines in the following order: Baselines east, -north,
            -longest.
        sta_labels : list of numpy.ndarray
        declination : float
        flag : bool
        symbol : str
        color : str
        sel_wl : float
        airmass_lim : float
        """
        # u, v = map(lambda x: x/sel_wl, uv_coords)
        u_coords, v_coords = uv_coord
        latitude_paranal = EarthLocation.of_site(
            "paranal").geodetic.lat.to(u.rad)
        hamax = np.arccos(abs((1./airmass_lim-np.sin(latitude_paranal)
                               * np.sin(declination))/(np.cos(latitude_paranal)
                                                       * np.cos(declination))))
        hour_angles = np.linspace(-hamax, hamax, 1000)
        u_coord_tracks, v_coord_tracks = self.calculate_uv_points(baselines,
                                                                  hour_angles,
                                                                  latitude_paranal,
                                                                  declination)
        # u, v = map(lambda x: x/sel_wl, ulvl)

        ax.plot(u_coord_tracks, v_coord_tracks, '-', color='grey', alpha=0.5)
        ax.plot(-u_coord_tracks, -v_coord_tracks, '-', color='grey', alpha=0.5)
        ax.plot([0.], [0.], '+k', markersize=5, markeredgewidth=2, alpha=0.5)

        ax.plot(u_coords, v_coords, symbol, color=color,
                markersize=10, markeredgewidth=3)
        ax.plot(-u_coords, -v_coords, symbol,
                color=color, markersize=10, markeredgewidth=3)
        ax.text(-u_coords-3.5, -v_coords-1.5, sta_label,
                fontsize="small", color='0', alpha=0.8)

    def plot_uv(self, ax: Axes, symbol: Optional[str] = "x",
                color: Optional[str | List[str]] = "b",
                sel_wl: Optional[float] = None,
                airmass_lim: Optional[float] = 2.) -> None:
        """Plots the (u, v)-coordinates and their corresponding tracks

        Parameters
        ----------
        ax : matplotlib.axes.Axes
        symbol : str, optional
            The symbol that markes the coordinates of the (u, v)-point.
        color : str or list of str, optional
            Set the color/colors of the uv-coords. In case of multiple
            (.fits)-files the colors can be specified as a list with entries
            for each file.
        sel_wl : float, optional
            The wavelength to convert meter to Mlambda. If None then no
            conversion is done.
        """
        for index, readout in enumerate(self.readouts):
            uv_coords = readout.oi_vis["UVCOORD"]
            flags = readout.oi_vis["FLAG"]
            sta_indices = readout.oi_vis["STA_INDEX"]
            sta_index = readout.oi_array["STA_INDEX"]
            sta_name = readout.oi_array["STA_NAME"]
            sta_xyz = readout.oi_array["STAXYZ"]

            baselines, sta_labels = [], []
            for uv_index, _ in enumerate(uv_coords):
                try:
                    baseline = sta_xyz[sta_indices[uv_index, 0] == sta_index][0]\
                        - sta_xyz[sta_indices[uv_index, 1] == sta_index][0]
                    sta_label = sta_name[sta_indices[uv_index, 0] == sta_index][0] + '-'\
                        + sta_name[sta_indices[uv_index, 1] == sta_index][0]
                except IndexError:
                    baseline, sta_label = [np.nan, np.nan, np.nan], ""

                baselines.append(baseline)
                sta_labels.append(sta_label)

            # TODO: Determine the maximum for the (u, v)-plot and then use that to set the
            # plot dimensions
            for uv_index, uv_coord in enumerate(uv_coords):
                self.make_uv_tracks(ax, uv_coord, baselines[uv_index],
                                    sta_labels[uv_index], readout.dec*np.pi/180,
                                    flags[uv_index], symbol,
                                    self.colors[index], sel_wl, airmass_lim)

            xlabel, ylabel = "$u$ (m)", "$v$ (m)"
            uv_max = np.nanmax(np.abs(uv_coords)) + 15

            ax.set_xlim((uv_max, -uv_max))
            ax.set_ylim((-uv_max, uv_max))
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

    def set_dataframe(self, labels: List[str], column: Column) -> DataFrame:
        """Prepares each row in a column as a pandas DataFrame"""
        return pd.DataFrame({label: array for label, array in zip(labels, column)})

    def make_component(self, data_name: str,
                       legend_format: Optional[str] = "long",
                       unwrap: Optional[bool] = False,
                       period: Optional[int] = 360) -> Union[Callable, DataFrame]:
        """Generates a pandas DataFrame that has all the plots' information

        Parameters
        ----------
        data_name: str
            The name of the data to be plotted. Determines the legend- and plot labels
        legend_format: str
            Sets the format of the legend: For all information set "long" and for only
            station names set "short"

        Returns
        -------
        component: list of either Callable or DataFrame
        """
        component = []
        for readout in self.readouts:
            if data_name == "flux":
                labels = readout.oi_array["TEL_NAME"]
                fluxes = readout.oi_flux["FLUXDATA"]
                if len(fluxes) == 1:
                    labels = ["Averaged"]
                sub_component = self.set_dataframe(labels, fluxes)
            elif data_name in ["vis", "vis2", "diff", "corrflux"]:
                station_names = readout.oi_vis["DELAY_LINE"]
                if legend_format == "long":
                    baselines = np.around(readout.oi_vis["BASELINE"], 2)
                    u_coords = readout.oi_vis["UVCOORD"][:, 0]
                    v_coords = readout.oi_vis["UVCOORD"][:, 1]
                    # TODO: Find out what this is exactly? Projected Baselines? Positional Angle?
                    pas = np.around(
                        (np.degrees(np.arctan2(v_coords, u_coords))-90)*-1, 2)
                    # TODO: Make the variables into mathrm
                    labels = [fr"{station_name} $B_p$={baseline} m $\phi={pa}^\circ$"
                              for station_name, baseline, pa in zip(station_names, baselines, pas)]
                else:
                    labels = station_names
                if data_name == "vis":
                    sub_component = self.set_dataframe(
                        labels, readout.oi_vis["VISAMP"])
                elif data_name == "diff":
                    # TODO: Make this into a function
                    diff_phases = readout.oi_vis["VISPHI"]
                    if unwrap:
                        diff_phases = np.unwrap(diff_phases, period=period)
                    sub_component = self.set_dataframe(labels, diff_phases)
                elif data_name == "vis2":
                    sub_component = self.set_dataframe(labels,
                                                       readout.oi_vis2["VIS2DATA"])
                else:
                    raise KeyError("No data-type of that data name exists!")
            elif data_name == "cphases":
                # TODO: Make this into a function
                cphases = readout.oi_t3["T3PHI"]
                if unwrap:
                    cphases = np.unwrap(cphases, period=period)
                sub_component = self.set_dataframe(readout.oi_t3["TRIANGLE"],
                                                   cphases)
            elif data_name == "uv":
                sub_component = self.plot_uv
            else:
                raise KeyError("Input data name cannot be queried!")

            if isinstance(sub_component, DataFrame):
                sub_component["lambda"] = readout.oi_wl["EFF_WAVE"].data.squeeze()

            component.append(sub_component)
        return component

    # TODO: Make unit support in this component
    def add_flux(self):
        """Adds the total flux(es) as a subplot"""
        self.components["Total Flux [Jy]"] = self.make_component("flux")
        return self

    # TODO: Add option to automatically set label depending on what the flux is
    def add_vis(self, corr_flux: Optional[bool] = False,
                legend_format: Optional[str] = "long"):
        """Adds the visibilities/correlated fluxes as a subplot"""
        label = "Correlated Flux [Jy]" if corr_flux else "Visibility"
        self.components[label] = self.make_component("vis", legend_format)
        return self

    def add_vis2(self, legend_format: Optional[str] = "long"):
        """Adds the squared visibilities as a subplot"""
        self.components["Squared Visibility"] = self.make_component(
            "vis2", legend_format)
        return self

    def add_cphases(self, unwrap: Optional[bool] = False, period: Optional[int] = 360):
        """Adds the closure phases as a subplot"""
        self.components[r"Closure phases [$^{\circ}$]"] =\
            self.make_component("cphases", unwrap=unwrap, period=period)
        return self

    def add_diff_phases(self, unwrap: Optional[bool] = False, period: Optional[int] = 360):
        """Adds the differential phases as a subplot"""
        self.components[r"Differential phases [$^{\circ}$]"] =\
            self.make_component("diff", unwrap=unwrap, period=period)
        return self

    def add_uv(self):
        """Adds the (u, v)-coordinates as a subplot"""
        self.components["$(u, v)$-coordinates"] = self.make_component("uv")
        return self

    def add_mosaic(self, unwrap: Optional[bool] = False,
                   legend_format: Optional[str] = "long"):
        """Combines multiple subplots to produce a mosaic plot"""
        self.add_uv().add_vis(False, legend_format).add_vis(True, legend_format)
        self.add_flux().add_cphases(unwrap).add_diff_phases(unwrap)
        return self

    def plot(self, save: Optional[bool] = False):
        """Combines the individual components into one plot.

        The size and dimension of the plot is automatically determined

        Parameters
        ----------
        save: bool, optional
            If toggled, saves the plot to the self.save_path file with the
            self.plot_name
        """
        columns = 1 if self.num_components == 1 else\
            (3 if self.num_components >= 3 else 2)
        rows = np.ceil(self.num_components/columns).astype(int)\
            if not self.num_components == 1 else 1
        to_px = 1/plt.rcParams["figure.dpi"]
        fig, axarr = plt.subplots(rows, columns,
                                  constrained_layout=True,
                                  figsize=(512*to_px*columns, 512*to_px*rows))
        # TODO: Add color and line support here
        if not self.num_components == 1:
            for ax, (name, component) in zip(axarr.flatten(), self.components.items()):
                for index, sub_component in enumerate(component):
                    if isinstance(sub_component, DataFrame):
                        sub_component.plot.line(x="lambda",
                                                xlabel=r"$\lambda$ [$\mathrm{\mu}$m]",
                                                ylabel=name, ax=ax, legend=True,
                                                linestyle=self.linestyles[index])
                        ax.legend(fontsize="x-small",
                                  loc="upper right", framealpha=0.5)
                    else:
                        sub_component(ax)
        else:
            name, component = map(
                lambda x: x[0], zip(*self.components.items()))
            for index, sub_component in enumerate(component):
                if isinstance(sub_component, DataFrame):
                    sub_component.plot.line(x="lambda",
                                            xlabel=r"$\lambda$ [$\mathrm{\mu}$m]",
                                            ylabel=name, ax=axarr, legend=True,
                                            linestyle=self.linestyles[index])
                    axarr.legend(fontsize="x-small",
                                 loc="upper right", framealpha=0.5)
                else:
                    sub_component(axarr)
        fig.tight_layout()

        if save:
            plt.savefig(str(self.save_path / self.plot_name), format="pdf")
        else:
            plt.show()
