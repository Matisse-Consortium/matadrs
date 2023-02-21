import pkg_resources
from pathlib import Path
from typing import Tuple, List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame
from astropy.table import Column

from ..utils.readout import ReadoutFits
from ..utils.tools import cprint


DATA_DIR = Path(pkg_resources.resource_filename("matadrs", "data"))

# Paranal lattitude in radians
LATITUDE = -24.62587 * np.pi / 180.

# TODO: Make this go away from class and just individual plotting functions? -> Maybe
# easier for paper relevant plots, but first make unified class that makes nice data
# reduction plots
class Plotter:
    """Class that plots models as well as reduced data

    Attributes
    ----------
    """

    def __init__(self, fits_file: Path,
                 flux_file: Optional[Path] = None,
                 plot_name: Optional[str] = None,
                 save_path: Optional[Path] = None) -> None:
        """Initialises the class instance"""
        self.readout = ReadoutFits(fits_file)
        self._band_mask = None

        # TODO: Make it that if multiple datasets are input that multiple plots are made
        # for the different bands -> See how to implement this

        if save_path is None:
            self.save_path = Path("").cwd()
        else:
            self.save_path = Path(save_path)

        # TODO: Improve this plot name giving. Only one fits file is respected in this
        if plot_name is None:
            self.plot_name = f"{Path(fits_file[0]).stem}.pdf"

        self.wl = self.readout.oi_wl["EFF_WAVE"].data[0]
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
                    self._band_mask = (wl > 1.6) & (wl < 1.8) & (wl > 3.) & (wl < 4.)
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

    def set_dataframe(self, labels: List[str], column: Column) -> DataFrame:
        """Prepares each row in a column as a pandas DataFrame"""
        return pd.DataFrame({label: array for label, array in zip(labels, column.data)})

    def make_component(self, data_name: str,
                       legend_format: Optional[str] = "long") -> DataFrame:
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
        df: DataFrame
        """
        if data_name == "flux":
            df = self.set_dataframe(self.readout.oi_flux["TEL_NAME"],
                                    self.readout.oi_flux["FLUXDATA"])
        elif (data_name == "vis") or (data_name == "vis2"):
            # TODO: Check if there are edge cases where oi_vis needs to be used
            station_names = self.readout.oi_vis["DELAY_LINE"]
            if legend_format == "long":
                baselines = np.around(self.readout.oi_vis["BASELINE"], 2)
                u_coords = self.readout.oi_vis["UVCOORD"][:, 0]
                v_coords = self.readout.oi_vis["UVCOORD"][:, 1]
                pas = np.around((np.degrees(np.arctan2(v_coords, u_coords))-90)*-1, 2)
                # TODO: Make the variables into mathrm
                labels = [fr"{station_name} $B_p$={baseline} m $\phi={pa}^\circ$"
                    for station_name, baseline, pa in zip(station_names, baselines, pas)]
            else:
                labels = station_names

            if data_name == "vis":
                df = self.set_dataframe(labels, self.readout.oi_vis["VISAMP"])
            else:
                df = self.set_dataframe(labels, self.readout.oi_vis2["VIS2DATA"])
            # TODO: Find out what this is exactly? Projected Baselines? Positional Angle?
        elif data_name == "cphase":
            df = self.set_dataframe(self.readout.oi_t3["TRIANGLE"],
                                    self.readout.oi_t3["T3PHI"])
        else:
            raise ValueError("Input data name cannot be queried!")
        df["lambda"] = self.wl
        # TODO: Why is this needed?
        return self.mask_dataframe(df, self.band_mask)

    def plot(self, save: Optional[bool] = False):
        """Makes the plot from the components

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
            # TODO: Make this more modular for future plots
            # TODO: Implement flux np.nan detection and not plotting it
            for ax, (name, dataframe) in zip(axarr.flatten(), self.components.items()):
                dataframe.plot(x="lambda", xlabel=r"$\lambda$ [$\mathrm{\mu}$m]",
                               ylabel=name, ax=ax, legend=True)
                plt.legend(fontsize=6)
        else:
            # TODO: Make this more modular for future plots
            name, dataframe = self.components[0].items()
            dataframe.plot(x="lambda", xlabel=r"$\lambda$ [$\mathrm{\mu}$m]",
                           ylabel=name, ax=axarr, legend=True)

        fig.tight_layout()

        if save:
            plt.savefig(str(self.save_path / self.plot_name), format="pdf")
        else:
            plt.show()
        plt.close()
        # ax.legend(loc=1, prop={'size': 6}, ncol=ncol)

    # TODO: Make somehow correlated flux and unit support in this component
    def add_flux(self):
        """Plots the flux """
        self.components["Flux [Jy]"] = self.make_component("flux")
        return self

    def add_vis(self, legend_format: Optional[str] = "long"):
        """Plots all the visibilities/correlated_fluxes in one plot """
        self.components["Vis"] = self.make_component("vis", legend_format)
        return self

    def add_vis2(self, legend_format: Optional[str] = "long"):
        """Plots all the visibilities squared in one plot"""
        self.components["Vis2"] = self.make_component("vis2", legend_format)
        return self

    def add_cphases(self):
        """Plots all the closure phases into one plot"""
        self.components["Closure phases [$^{\circ}$]"] = self.make_component("cphase")
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

    # TODO: Rename function at a future point
    def calculate_uv_points(self, baselines: List[float],
                            hour_angle: np.ndarray[float]) -> Tuple[np.ndarray]:
        """Calculates the earth rotation (synthesis) for the uv-point corresponding to the
        baselines for the input hour angle(s)

        Parameters
        -----------
        baselines: List[np.ndarray]
        hour_angle: np.ndarray[float]

        Returns
        -------
        u_coords: np.ndarray
        v_coords: np.ndarray
        """
        baseline_east, baseline_north, baseline_longest = baselines
        dec_rad = self.readout.dec*np.pi/180

        u_coords = baseline_east * np.cos(hour_angle) - baseline_north * np.sin(LATITUDE)\
            * np.sin(hour_angle) + baseline_longest * np.cos(LATITUDE) * np.sin(hour_angle)
        v_coords = baseline_east * np.sin(dec_rad) * np.sin(hour_angle)\
            + baseline_north * (np.sin(LATITUDE) * np.sin(dec_rad) * np.cos(hour_angle)\
                + np.cos(LATITUDE) * np.cos(dec_rad)) - baseline_longest * \
                (np.cos(LATITUDE) * np.sin(dec_rad) * np.cos(hour_angle)\
                    - np.sin(LATITUDE) * np.cos(dec_rad))
        return u_coords, v_coords

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
            The baselines in the following order: Baselines east, -north, -longest
        sta_labels: List[np.ndarray]
        flag: bool
        symbol: str
        color: str
        """
        dec_rad = self.readout.dec*np.pi/180
        # u, v = map(lambda x: x/sel_wl, uv_coords)
        u_coords, v_coords = uv_coord

        # if not color:
        #     if np.all(flag) == 'True':
        #         color = 'r'
        #     else:
        #         color = 'g'

        hamax = np.arccos(abs((1. / airmass_lim - np.sin(LATITUDE) * np.sin(dec_rad))/\
                              (np.cos(LATITUDE) * np.cos(dec_rad))))
        hour_angles = np.linspace(-hamax, hamax, 1000)
        u_coord_tracks, v_coord_tracks = self.calculate_uv_points(baselines, hour_angles)
        # u, v = map(lambda x: x/sel_wl, ulvl)

        ax.plot(u_coord_tracks, v_coord_tracks, '-', color='grey', alpha=0.5)
        ax.plot(-u_coord_tracks, -v_coord_tracks, '-', color='grey', alpha=0.5)
        ax.plot([0.], [0.], '+k', markersize=5, markeredgewidth=2, alpha=0.5)

        ax.plot(u_coords, v_coords, symbol, color=color,
                markersize=10, markeredgewidth=3)
        ax.plot(-u_coords, -v_coords, symbol,
                color=color, markersize=10, markeredgewidth=3)
        ax.text(-u_coords-7, -v_coords-3, sta_label, color='0', alpha=0.8)

    def plot_uv(self, ax, symbol: Optional[str] = "x",
                color: Optional[str] = "b",
                sel_wl: Optional[float] = None,
                airmass_lim: Optional[float] = 2.):
        """Makes the (u, v)-plots

        Parameters
        ----------
        ax
        symbol: str, optional
        color: str, optional
        sel_wl: float, optional
        """
        uv_coords = self.readout.oi_vis2["UVCOORD"]
        flags = self.readout.oi_vis2["FLAG"]
        sta_indices = self.readout.oi_vis2["STA_INDEX"]
        sta_index = self.readout.oi_array["STA_INDEX"]
        sta_name = self.readout.oi_array["STA_NAME"]
        sta_xyz = self.readout.oi_array["STAXYZ"]

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

        # umax = np.nanmax(np.abs(u))
        # vmax = np.nanmax(np.abs(v))
        # if self.readout.mindexd is None:
        #     mindexd = np.amin(self.readout.oi_vis2["MindexD"][0])
        # else:
        # mindexd = self.readout.mindexd
        # rel_time = (self.readout.vis2["MindexD"] - mindexd) * 24.0 * 3600.0  # (s)
        # dic['TREL'] = rel_time[0]

        for index, uv_coord in enumerate(uv_coords):
            self.make_uv_tracks(ax, uv_coord, baselines[index],
                                sta_labels[index], flags[index],
                                symbol, color, sel_wl, airmass_lim)

        xlabel = "$u$ (m)"
        ylabel = "$v$ (m)"
        # xlabel = r"$u$ ($M\lambda$)"
        # ylabel = r"$v$ ($M\lambda$)"
        ax.set_xlim((150, -150))
        ax.set_ylim((-150, 150))
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


if __name__ == ('__main__'):
    fits_files = ["HD_163296_2019-03-23T08_41_19_N_TARGET_FINALCAL_INT.fits",
                  "HD_163296_2019-03-23T08_41_19_L_TARGET_FINALCAL_INT.fits",
                  "HD_163296_2019-05-06T08_19_51_L_TARGET_FINALCAL_INT.fits"]
    fits_files = [DATA_DIR / "tests" / fits_file for fits_file in fits_files]
    plot_fits = Plotter(fits_files[1])
    plot_fits.add_cphases().add_vis("short").plot()
