from pathlib import Path
from typing import Callable, Tuple, List, Union, Optional

import astropy.units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame
from astropy.coordinates import EarthLocation
from astropy.table import Column

from ..utils.data_prep import DataPrep


# Paranal lattitude in radians
LATITUDE = -24.62587 * np.pi / 180.


# TODO: Make this go away from class and just individual plotting functions? -> Maybe
# easier for paper relevant plots, but first make unified class that makes nice data
# reduction plots
class Plotter:
    """Class that plots models as well as reduced data

    Parameters
    ----------

    Attributes
    ----------

    Methods
    -------
    """

    def __init__(self, fits_files: Path | List[Path],
                 flux_file: Optional[Path] = None,
                 plot_name: Optional[str] = None,
                 save_path: Optional[Path] = None) -> None:
        """Initialises the class instance"""
        self.data_prep = DataPrep(fits_files)
        self._band_mask = None

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

        self.wl = self.data_prep.oi_wl["EFF_WAVE"].data[0]
        self.components = {}

    @property
    def number_of_plots(self):
        return len(self.components)

    # TODO: Make this more modular for multiple files -> Right now only one wl is
    # respected
    @property
    def band_mask(self):
        """The masking for the bands to make the plots visually readable"""
        if self._band_mask is None:
            wl = pd.DataFrame({"wl": self.wl})
            if np.any(wl < 7.):
                if np.any(wl < 2.) and np.any(wl > 3.):
                    self._band_mask = (wl > 1.6) & (
                        wl < 1.8) & (wl > 3.) & (wl < 4.)
                elif np.any(wl < 2.):
                    self._band_mask = (wl > 1.6) & (wl < 1.8)
                else:
                    self._band_mask = (wl > 3.) & (wl < 4.)
            else:
                self._band_mask = (wl > 8.5) & (wl < 12.5)
        return ~self._band_mask["wl"]

    # HACK: This should not need to be a thing, but it doesn't work otherwise
    def mask_dataframe(self, df: DataFrame, mask: np.ndarray) -> None:
        """Iterates through each row to mask a DataFrame"""
        for column in df.columns:
            df[column] = df[column].mask(mask)
        return df

    # TODO: Rename function at a future point
    def calculate_uv_points(self, baselines: List[float],
                            hour_angle: np.ndarray[float],
                            latitude: u.rad) -> Tuple[np.ndarray]:
        """Calculates the earth rotation (synthesis) for the uv-point
        corresponding to the baselines for the input hour angle(s)

        Parameters
        -----------
        baselines: List[float]
            The baselines in the following order: Baselines east, -north, -longest
        hour_angle: np.ndarray[float]
        latitude: u.rad
            The latitude of the site

        Returns
        -------
        u_coords: np.ndarray
        v_coords: np.ndarray
        """
        baseline_east, baseline_north, baseline_longest = baselines
        dec_rad = self.data_prep.dec*np.pi/180

        u_coords = baseline_east * np.cos(hour_angle) - baseline_north * np.sin(latitude)\
            * np.sin(hour_angle) + baseline_longest * np.cos(latitude) * np.sin(hour_angle)
        v_coords = baseline_east * np.sin(dec_rad) * np.sin(hour_angle)\
            + baseline_north * (np.sin(latitude) * np.sin(dec_rad) * np.cos(hour_angle)
                                + np.cos(latitude) * np.cos(dec_rad)) - baseline_longest * \
            (np.cos(latitude) * np.sin(dec_rad) * np.cos(hour_angle)
             - np.sin(latitude) * np.cos(dec_rad))
        return u_coords, v_coords

    # TODO: Add proper docs
    def make_uv_tracks(self, ax, uv_coord: np.ndarray[float],
                       baselines: List[np.ndarray],
                       sta_label: List[np.ndarray], flag: bool,
                       symbol: str, color: str, sel_wl: float,
                       airmass_lim: float) -> None:
        """This function was written by Jozsef Varga (from menEWS: menEWS_plot.py).

        From coordinate + ha (range), calculate uv tracks

        Parameters
        ----------
        uv_coords: np.ndarray[float]
        baselines: List[np.ndarray]
            The baselines in the following order: Baselines east, -north,
            -longest
        sta_labels: List[np.ndarray]
        flag: bool
        symbol: str
        color: str
        sel_wl: float
        airmass_lim: float
        """
        # TODO: Find a way how to make this calculation possible for different datasets
        dec_rad = self.data_prep.dec*np.pi/180
        # u, v = map(lambda x: x/sel_wl, uv_coords)
        u_coords, v_coords = uv_coord
        latitude_paranal = EarthLocation.of_site(
            "paranal").geodetic.lat.to(u.rad)
        hamax = np.arccos(abs((1./airmass_lim-np.sin(latitude_paranal)
                               * np.sin(dec_rad))/(np.cos(latitude_paranal)
                                                   * np.cos(dec_rad))))
        hour_angles = np.linspace(-hamax, hamax, 1000)
        u_coord_tracks, v_coord_tracks = self.calculate_uv_points(baselines,
                                                                  hour_angles,
                                                                  latitude_paranal)
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

    def plot_uv(self, ax, symbol: Optional[str] = "x",
                color: Optional[str | List[str]] = "b",
                sel_wl: Optional[float] = None,
                airmass_lim: Optional[float] = 2.) -> None:
        """Plots the (u, v)-coordinates and their corresponding tracks

        Parameters
        ----------
        ax
        symbol: str, optional
            The symbol that markes the coordinates of the (u, v)-point
        color: str | List[str], optional
            Set the color/colors of the uv-coords. In case of multiple
            (.fits)-files the colors can be specified as a list with entries
            for each file
        sel_wl: float, optional
            The wavelength to convert meter to Mlambda. If None then no
            conversion is done
        """
        uv_coords = self.data_prep.oi_vis["UVCOORD"]
        flags = self.data_prep.oi_vis["FLAG"]
        sta_indices = self.data_prep.oi_vis["STA_INDEX"]
        sta_index = self.data_prep.oi_array["STA_INDEX"]
        sta_name = self.data_prep.oi_array["STA_NAME"]
        sta_xyz = self.data_prep.oi_array["STAXYZ"]

        baselines, sta_labels = [], []
        for index, _ in enumerate(uv_coords):
            try:
                baseline = sta_xyz[sta_indices[index, 0] == sta_index][0]\
                    - sta_xyz[sta_indices[index, 1] == sta_index][0]
                sta_label = sta_name[sta_indices[index, 0] == sta_index][0] + '-'\
                    + sta_name[sta_indices[index, 1] == sta_index][0]
            except IndexError:
                baseline, sta_label = [np.nan, np.nan, np.nan], ""

            baselines.append(baseline)
            sta_labels.append(sta_label)

        # TODO: Determine the maximum for the (u, v)-plot and then use that to set the
        # plot dimensions

        # if self.readout.mindexd is None:
        #     mindexd = np.amin(self.readout.oi_vis2["MindexD"][0])
        # else:
        # mindexd = self.readout.mindexd
        # rel_time = (self.readout.vis2["MindexD"] - mindexd) * 24.0 * 3600.0  # (s)
        # dic['TREL'] = rel_time[0]

        fits_index = 0
        for index, uv_coord in enumerate(uv_coords):
            if index+1 % 6 == 0:
                fits_index += 1
            self.make_uv_tracks(ax, uv_coord, baselines[index],
                                sta_labels[index], flags[index],
                                symbol, color, sel_wl, airmass_lim)

        xlabel, ylabel = "$u$ (m)", "$v$ (m)"
        # xlabel, ylabel = r"$u$ ($M\lambda$)", r"$v$ ($M\lambda$)"
        uv_max = np.nanmax(np.abs(uv_coords)) + 15

        ax.set_xlim((uv_max, -uv_max))
        ax.set_ylim((-uv_max, uv_max))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        # plotmax = 1.3*np.amax([umax, vmax])

        # plot_title = dic['TARGET'] + "\n" + "date: " + dic['DATE-OBS'] + "\n" + "TPL start: " + dic['TPL_START'] + "\n" + dic['CATEGORY'] + ' ' +\
        #     dic['BAND'] + ' ' + dic['DISPNAME'] #+ ' ' + dic['BCD1'] + '-' + dic['BCD2']
        # if math.isnan(B_lim[0]):
        #     xlim = (+plotmax / sel_wl, -plotmax / sel_wl)
        #     ylim = (-plotmax / sel_wl, +plotmax / sel_wl)
        # else:
        #     xlim = (+B_lim[1] / sel_wl, -B_lim[1] / sel_wl)
        #     ylim = (-B_lim[1] / sel_wl, +B_lim[1] / sel_wl)
        # if plot_Mlambda == True:
        # plot_config(xlabel, ylabel,plot_title, ax, dic,
        #             ylim=ylim,xlim=xlim,plot_legend=False,annotate=annotate)

    def set_dataframe(self, labels: List[str], column: Column) -> DataFrame:
        """Prepares each row in a column as a pandas DataFrame"""
        if isinstance(column, Column):
            data = column.data
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
        component: Callable | DataFrame
        """
        if data_name == "flux":
            labels = self.data_prep.oi_array["TEL_NAME"]
            flux = self.data_prep.oi_flux["FLUXDATA"]
            if len(flux) == 1:
                labels = ["Averaged"]
            component = self.set_dataframe(labels, flux)
        elif data_name in ["vis", "vis2", "diff", "corrflux"]:
            station_names = self.data_prep.oi_vis["DELAY_LINE"]
            if legend_format == "long":
                baselines = np.around(self.data_prep.oi_vis["BASELINE"], 2)
                u_coords = self.data_prep.oi_vis["UVCOORD"][:, 0]
                v_coords = self.data_prep.oi_vis["UVCOORD"][:, 1]
                # TODO: Find out what this is exactly? Projected Baselines? Positional Angle?
                pas = np.around(
                    (np.degrees(np.arctan2(v_coords, u_coords))-90)*-1, 2)
                # TODO: Make the variables into mathrm
                labels = [fr"{station_name} $B_p$={baseline} m $\phi={pa}^\circ$"
                          for station_name, baseline, pa in zip(station_names, baselines, pas)]
            else:
                labels = station_names
            if data_name == "vis":
                component = self.set_dataframe(
                    labels, self.data_prep.oi_vis["VISAMP"])
            elif data_name == "diff":
                # TODO: Make this into a function
                diff_phases = self.data_prep.oi_vis["VISPHI"]
                if unwrap:
                    diff_phases = np.unwrap(diff_phases, period=period)
                component = self.set_dataframe(labels, diff_phases)
            elif data_name == "corrflux":
                try:
                    component = self.set_dataframe(
                        labels, self.data_prep.oi_cfx["VISAMP"])
                except KeyError:
                    return self.make_component("vis", legend_format)
            elif data_name == "vis2":
                component = self.set_dataframe(
                    labels, self.data_prep.oi_vis2["VIS2DATA"])
            else:
                raise KeyError("No data-type of that data name exists!")
        elif data_name == "cphases":
            # TODO: Make this into a function
            cphases = self.data_prep.oi_t3["T3PHI"]
            if unwrap:
                cphases = np.unwrap(cphases, period=period)
            component = self.set_dataframe(
                self.data_prep.oi_t3["TRIANGLE"], cphases)
        elif data_name == "uv":
            component = self.plot_uv
        else:
            raise KeyError("Input data name cannot be queried!")
        if isinstance(component, DataFrame):
            # HACK: This is needed right now, but it should really be doable with easier
            # operations
            component["lambda"] = self.wl
            return self.mask_dataframe(component, self.band_mask)
        return component

    def plot(self, save: Optional[bool] = False):
        """Makes the plot from the omponents

        Parameters
        ----------
        save: bool, optional
            If toggled, saves the plot to the self.save_path file with the
            self.plot_name
        """
        # TODO: Handle save structure differently -> with save_path and so
        columns = 1 if self.number_of_plots == 1 else\
            (3 if self.number_of_plots >= 3 else 2)
        rows = np.ceil(self.number_of_plots/columns).astype(int)\
            if not self.number_of_plots == 1 else 1
        to_px = 1/plt.rcParams["figure.dpi"]
        fig, axarr = plt.subplots(rows, columns,
                                  constrained_layout=True,
                                  figsize=(512*to_px*columns, 512*to_px*rows))
        if not self.number_of_plots == 1:
            # TODO: Implement if flux is np.nan and not plotting it in that case
            for ax, (name, component) in zip(axarr.flatten(), self.components.items()):
                if isinstance(component, DataFrame):
                    component.plot(x="lambda", xlabel=r"$\lambda$ [$\mathrm{\mu}$m]",
                                   ylabel=name, ax=ax, legend=True)
                    ax.legend(fontsize="x-small",
                              loc="upper right", framealpha=0.5)
                else:
                    component(ax)
        else:
            name, component = zip(*self.components.items())
            if isinstance(component, DataFrame):
                component.plot(x="lambda", xlabel=r"$\lambda$ [$\mathrm{\mu}$m]",
                               ylabel=name, ax=axarr, legend=True)
                axarr.legend(fontsize="x-small",
                             loc="upper right", framealpha=0.5)
            else:
                component[0](axarr)

        fig.tight_layout()

        if save:
            plt.savefig(str(self.save_path / self.plot_name), format="pdf")
        else:
            plt.show()

    # TODO: Make somehow correlated flux and unit support in this component
    def add_flux(self):
        """Plots the flux """
        self.components["Total Flux [Jy]"] = self.make_component("flux")
        return self

    # TODO: Add option to remove the correlated flux if the component is empty
    def add_vis(self, corr_flux: Optional[bool] = False,
                legend_format: Optional[str] = "long"):
        """Plots all the visibilities fluxes in one plot """
        if corr_flux:
            self.components["Correlated Flux [Jy]"] = self.make_component("corrflux",
                                                                          legend_format)
        else:
            self.components["Visibility"] = self.make_component(
                "vis", legend_format)
        return self

    def add_vis2(self, legend_format: Optional[str] = "long"):
        """Plots all the visibilities squared in one plot"""
        self.components["Squared Visibility"] = self.make_component(
            "vis2", legend_format)
        return self

    def add_diff_phases(self, unwrap: Optional[bool] = False, period: Optional[int] = 360):
        """Plots all the differential phases into one plot"""
        self.components["Differential phases [$^{\circ}$]"] =\
            self.make_component("diff", unwrap=unwrap, period=period)
        return self

    def add_cphases(self, unwrap: Optional[bool] = False, period: Optional[int] = 360):
        """Plots all the closure phases into one plot"""
        self.components["Closure phases [$^{\circ}$]"] =\
            self.make_component("cphases", unwrap=unwrap, period=period)
        return self

    def add_uv(self):
        """Plots the (u, v)-coordinates"""
        self.components["$(u, v)$-coordinates"] = self.make_component("uv")
        return self

    def add_mosaic(self, unwrap: Optional[bool] = False,
                   legend_format: Optional[str] = "long"):
        """Prepares a combined mosaic-plot"""
        self.add_uv().add_vis(False, legend_format).add_vis(True, legend_format)
        self.add_flux().add_cphases(unwrap).add_diff_phases(unwrap)
        return self

    # def get_plot_linestyle(self, already_chosen_linestyles: List):
        # """Gets a linestyle, which is different from the already chosen one

        # Parameters
        # ----------
        # already_chosen_linestyles: List
        # A list that contains the linestyles that have already been applied
        # """
        # linestyles = ["-", "--", "-.", ":"]
        # for linestyle in linestyles:
        # if linestyle not in already_chosen_linestyles:
        # return linestyle

        # if all(linestyle in linestyles in already_chosen_linestyles):
        # raise RuntimeError("No linestyle to pick as all have been picked already!")

    # def plot_vis24baseline(self, ax, do_fit: Optional[bool] = True) -> None:
        # """ Plot the mean visibility for one certain wavelength and fit it with a
        # gaussian and airy disk

        # Parameters
        # ----------
        # ax
        # The matplotlib axis to be used for the plot
        # do_fit: bool, optional
        # Switches the fit mechanic on or off. DEFAULT "True"
        # """
        # ax.plot(self.baseline_distances, self.mean_bin_vis2, ls='None', marker='o')

        # # TODO: Implement fit here
        # if do_fit:
        # ...

        # ax.set_xlabel(fr'uv-distance [m] at $\lambda_0$={self.mean_wl:.2f} $\mu m$')
        # ax.set_ylabel(r'$\bar{V}$')

    # def plot_waterfall(self, ax) -> None:
        # """Plots a waterfall with the mean wavelength for the different baselines"""
        # for i in range(6):
        # ax.errorbar(self.wl[self.si:self.ei]*1e06,
        # self.vis2data[i][self.si:self.ei],
        # yerr=np.nanstd(self.vis2data[i][self.si:self.ei]),
        # label=self.tel_vis2[i], ls='None', fmt='o')
        # ax.set_xlabel(r'wl [micron]')
        # ax.set_ylabel('vis2')
        # ax.legend(loc='best')
