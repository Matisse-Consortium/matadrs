import os
from pathlib import Path
from warnings import warn
from typing import List, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame
from astropy.table import Table, Column

# TODO: Find way to make this into a complete module -> More pythonic!
from data_prep import DataPrep

# TODO: Make either model or fourier transform carry more info like the name of
# the plot or similar -> Work more with classes
# TODO: Facilitate the plotting by getting the data from the DataPrep class and putting it
# into a pandas frame to easily plot it

DATA_DIR = Path(__file__).parent.parent.parent / "data"


class Plotter:
    """Class that plots models as well as vis-, t3phi- and uv-data

    Attributes
    ----------
    """
    def __init__(self, fits_files: List[Path],
                 flux_files: Optional[List[Path]] = None,
                 plot_name: Optional[str] = None,
                 save_path: Optional[Path] = None) -> None:
        """Initialises the class instance"""
        # TODO: Implement L- or N-band check. Fix this
        self.lband = False
        # TODO: Make it that if multiple datasets are input that multiple plots are made
        # for the different bands -> See how to implement this
        if len(fits_files) > 1:
            warn("There is more than one (.fits)-file detected."\
                 " WARNING: Only data of SAME BAND can be unified into one plot!")
        self.data_prep = DataPrep(fits_files)

        # TODO: Improve this plot name giving. Only one fits file is respected in this
        if plot_name is None:
            self.plot_name = "combined_files.pdf" if len(fits_files) > 1\
                else f"{fits_files[0].name.split('.')[0]}.pdf"
        else:
            self.plot_name = plot_name

        self.save_path = save_path
        self.components = {}

    @property
    def number_of_plots(self):
        return len(self.components)

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

    def set_dataframe(self, labels: List[str], column: Column) -> DataFrame:
        """Prepares each row in a column as a pandas DataFrame"""
        return pd.DataFrame({label: array for label, array in zip(labels, column.data)})

    # TODO: Settings for dataframe cuts if lband and general cuts. IMPLEMENT This!
            # if 'MATISSE' in dic['INSTRUMENT']:
                # if math.isnan(wl_lim[0]) or math.isnan(wl_lim[1]):
                    # M_idx = np.logical_not(np.logical_and(x > 1.6,x < 1.8)+np.logical_and(x > 3.0,x < 4.0)+np.logical_and(x > 8.0,x < 12.5) )
                # else:
                    # M_idx = np.logical_not(np.logical_and(x > 1.6,x < 1.8)+np.logical_and(x > 3.0,x < 4.0)+np.logical_and(x > 8.0,x < 12.5)+\
                    # np.logical_and(x > wl_lim[0],x < wl_lim[1]) )

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
            df = self.set_dataframe(self.data_prep.oi_flux["TEL_NAME"],
                                    self.data_prep.oi_flux["FLUXDATA"])
        elif (data_name == "vis") or (data_name == "vis2"):
            # TODO: Check if there are edge cases where oi_vis needs to be used
            station_names = self.data_prep.oi_vis["DELAY_LINE"]
            if legend_format == "long":
                baselines = np.around(self.data_prep.oi_vis["BASELINE"], 2)
                u_coords = self.data_prep.oi_vis["UVCOORD"][:, 0]
                v_coords = self.data_prep.oi_vis["UVCOORD"][:, 1]
                pas = np.around((np.degrees(np.arctan2(v_coords, u_coords))-90)*-1, 2)
                # TODO: Make the variables into mathrm
                labels = [fr"{station_name} $B_p$={baseline} m $\phi={pa}^\circ$"\
                        for station_name, baseline, pa in zip(station_names, baselines, pas)]
            else:
                labels = station_names

            if data_name == "vis":
                df = self.set_dataframe(labels, self.data_prep.oi_vis["VISAMP"])
            else:
                df = self.set_dataframe(labels, self.data_prep.oi_vis2["VIS2DATA"])
            # TODO: Find out what this is exactly? Projected Baselines? Positional Angle?
        elif data_name == "cphase":
            df =  self.set_dataframe(self.data_prep.oi_t3["TRIANGLE"],
                                     self.data_prep.oi_t3["T3PHI"])
        else:
            raise ValueError("Input data name cannot be queried!")
        # TODO: Make this more modular for multiple files
        df["lambda"] = self.data_prep.oi_wl["EFF_WAVE"].data[0]
        return df

    def plot(self, save: Optional[bool] = False):
        """Makes the plot from the components

        Parameters
        ----------
        save: bool, optional
            If toggled, saves the plot to the self.save_path file with the
            self.plot_name
        """
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
                               ylabel=name, ax=ax)
        else:
            # TODO: Make this more modular for future plots
            name, dataframe = self.components[0].items()
            dataframe.plot(x="lambda", xlabel=r"$\lambda$ [$\mathrm{\mu}$m]",
                           ylabel=name, ax=axarr)

        fig.tight_layout()

        if save:
            plt.savefig(self.save_path / self.plot_name, format="pdf")
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

    def add_cphase(self):
        """Plots all the closure phases into one plot"""
        self.components["Closure phases [$^{\circ}$]"] = self.make_component("cphase")
        return self

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

    # def plot_uv(self, ax) -> None:
        # """Plots the uv-coordinates with an orientational compass and the synthethis

        # Parameters
        # ----------
        # ax
            # The axis anchor of matplotlib.pyplot

        # Returns
        # -------
        # None
        # """
        # ax.scatter(self.ucoords, self.vcoords)
        # ax.scatter(-self.ucoords, -self.vcoords)
        # ax.set_xlim([150, -150])
        # ax.set_ylim([-150, 150])
        # ax.set_ylabel('v [m]')
        # ax.set_xlabel('u [m]')

        # # Compass for the directions
        # cardinal_vectors = [(0,1), (0,-1), (1,0), (-1,0)]
        # cardinal_colors  = ['black', 'green', 'blue', 'red']
        # cardinal_directions = ['N', 'S', 'W', 'E']
        # arrow_size, head_size = 40, 10
        # x, y = (-85, 85)

        # for vector, color, direction in zip(cardinal_vectors,
                                            # cardinal_colors, cardinal_directions):
            # dx, dy = vector[0]*arrow_size, vector[1]*arrow_size
            # if vector[0] == 0:
                # ax.text(x-dx-5, y+dy, direction)
            # if vector[1] == 0:
                # ax.text(x-dx, y+dy+5, direction)
            # arrow_args = {"length_includes_head": True,
                          # "head_width": head_size, "head_length": head_size,
                          # "width": 1, "fc": color, "ec": color}
            # ax.arrow(x, y, dx, dy, **arrow_args)


# def rotation_synthesis_uv(inp):
    # """This function was written by Jozsef Varga (from menEWS: menEWS_plot.py).

    # Calculates uv-point corresponding to inp (see "get_header_info"),
    # for hour angle(s) (ha)
    # """
    # ra, dec, BE, BN, BL, base = inp
    # paranal_lat = -24.62587 * np.pi / 180.

    # u = BE * np.cos(ha) -\
            # BN * np.sin(lat) * np.sin(ha) + BL * np.cos(lat) * np.sin(ha)
    # v = BE * np.sin(dec) * np.sin(ha) +\
            # BN * (np.sin(lat) * np.sin(dec) * np.cos(ha) +\
                  # np.cos(lat) * np.cos(dec)) - BL * \
        # (np.cos(lat) * np.sin(dec) * np.cos(ha)- np.sin(lat) * np.cos(dec))
    # return u, v


# def make_uv_tracks(uv, inp, flag, ax, bases=[], symbol='x',color='',
    # print_station_names=True,sel_wl=1.0,plot_Mlambda=False):
    # """This function was written by Jozsef Varga (from menEWS: menEWS_plot.py).

    # From coordinate + ha (range), calculate uv tracks"""

    # ra, dec, BE, BN, BL, base = inp
    # paranal_lat = -24.62587 * np.pi / 180.
    # mlim = 2.0  # airmass limit for tracks

    # if plot_Mlambda == True:
        # u, v = map(lambda x: x/sel_wl, uv)
    # else:
        # u, v = uv

    # if not color:
        # if np.all(flag) == 'True':
            # color = 'r'
        # else:
            # color = 'g'

    # if base not in bases:
        # hamax = np.arccos(abs((1. / mlim - np.sin(lat) * np.sin(dec)) / \
                              # (np.cos(lat) * np.cos(dec))))
        # harng = np.linspace(-hamax, hamax, 1000)

        # ul, vl = ulvl = calculate_uv_points(inp, harng)
        # if plot_Mlambda == True:
            # u, v = map(lambda x: x/sel_wl, ulvl)

        # ax.plot(ul, vl, '-', color='grey',alpha=0.5)
        # ax.plot(-ul, -vl, '-', color='grey',alpha=0.5)
        # ax.plot([0.], [0.], '+k', markersize=5, markeredgewidth=2,alpha=0.5)

        # if print_station_names:
            # ax.text(-u-7, -v-3, base, color='0',alpha=0.8)
        # bases.append(base)

    # ax.plot(u, v, symbol, color=color, markersize=10, markeredgewidth=3)
    # ax.plot(-u, -v, symbol, color=color, markersize=10, markeredgewidth=3)

    # return bases


# def make_uv_plot(dic,ax,verbose=False,annotate=True,B_lim=(np.nan,np.nan),figsize=(5,5),
    # color='',print_station_names=True,sel_wl=1.0,plot_Mlambda=False):
    # """This function was written by Jozsef Varga (from menEWS: menEWS_plot.py)"""
    # if plot_Mlambda==False:
        # sel_wl = 1.0
    # try:
        # u = dic['VIS2']['U']
        # v = dic['VIS2']['V']
        # flag = dic['VIS2']['FLAG']
        # sta_index = dic['VIS2']['STA_INDEX']
        # mjd = dic['VIS2']['MJD']
    # except KeyError as e:
        # if verbose: print(e)
        # u = [0.0]
        # v = [0.0]
        # flags = [False]
        # sta_index = []
        # mjd = [0.0]

    # uvs = []
    # inps = []
    # flags = []
    # umax = []
    # vmax = []
    # for j in range(len(u)):
        # uvs.append([u[j],v[j]])
        # try:
            # BE, BN, BL = dic['STAXYZ'][sta_index[j, 0] == dic['STA_INDEX']][0] - \
                # dic['STAXYZ'][sta_index[j, 1] == dic['STA_INDEX']][0]
            # sta_label= dic['STA_NAME'][sta_index[j, 0] == dic['STA_INDEX']][0] + '-' + \
                        # dic['STA_NAME'][sta_index[j, 1] == dic['STA_INDEX']][0]
        # except IndexError as e:
            # print('make_uv_plot STA_INDEX error.')
            # print(e)
            # BE, BN, BL = [np.nan,np.nan,np.nan]
            # sta_label= ''
        # inps.append( [dic['RA'] * np.pi / 180., dic['DEC'] * np.pi / 180., BE, BN, BL, sta_label]  )
        # flags.append(flag[j])
    # bases = []
    # umax = np.nanmax(np.abs(u))
    # vmax = np.nanmax(np.abs(v))
    # if not (dic['MJD-OBS']):
        # dic['MJD-OBS'] = np.amin(mjd[0])
    # try:
        # rel_time = (mjd - dic['MJD-OBS']) * 24.0 * 3600.0  # (s)
        # dic['TREL'] = rel_time[0]

        # for k, uv in enumerate(uvs):
            # bases = make_uv_tracks(uv, inps[k], flags[k],ax, bases,
            # color=color,print_station_names=print_station_names,
            # sel_wl=sel_wl,plot_Mlambda=plot_Mlambda)

        # if plot_Mlambda == False:
            # xlabel ='$u$ (m)'
            # ylabel ='$v$ (m)'
        # else:
            # xlabel ='$u$ ($M\lambda$)'
            # ylabel ='$v$ ($M\lambda$)'
        # ax.set_xlim((130, -130))
        # ax.set_ylim((-130, 130))
        # plotmax = 1.3*np.amax([umax,vmax])

        # plot_title = dic['TARGET'] + "\n" + "date: " + dic['DATE-OBS'] + "\n" + "TPL start: " + dic['TPL_START'] + "\n" + dic['CATEGORY'] + ' ' +\
            # dic['BAND'] + ' ' + dic['DISPNAME'] #+ ' ' + dic['BCD1'] + '-' + dic['BCD2']
        # if math.isnan(B_lim[0]):
            # xlim = (+plotmax/ sel_wl,-plotmax/ sel_wl)
            # ylim = (-plotmax/ sel_wl,+plotmax/ sel_wl)
        # else:
            # xlim = (+B_lim[1]/ sel_wl,-B_lim[1]/ sel_wl)
            # ylim = (-B_lim[1]/ sel_wl,+B_lim[1]/ sel_wl)
        # #if plot_Mlambda == True:
        # plot_config(xlabel, ylabel,plot_title, ax, dic,
                    # ylim=ylim,xlim=xlim,plot_legend=False,annotate=annotate)
    # except TypeError as e:
        # if verbose: print('Unable to plot ' + 'uv')
        # if verbose: print(e)
        # return 1

    # return 0


if __name__ == ('__main__'):
    fits_files = ["HD_163296_2019-03-23T08_41_19_N_TARGET_FINALCAL_INT.fits"]
                  # "HD_163296_2019-03-23T08_41_19_L_TARGET_FINALCAL_INT.fits",
                  # "HD_163296_2019-05-06T08_19_51_L_TARGET_FINALCAL_INT.fits"]
    fits_files = [DATA_DIR / "tests" / fits_file for fits_file in fits_files]
    plot_fits = Plotter(fits_files)
    plot_fits.add_cphase().add_vis().plot()


