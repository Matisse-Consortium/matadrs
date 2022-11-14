import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from glob import glob
from pathlib import Path
from warnings import warn
from astropy.io import fits
from typing import Any, Dict, List, Optional, Union

from do_readout import ReadoutFits

# TODO: Make either model or fourier transform carry more info like the name of
# the plot or similar -> Work more with classes
# TODO: Remove error bars from plots
# TODO: Rework dump function into pickle dump

# TODO: Implement plot save func

# TODO: Make standalone method that makes plot plot immediately
# TODO: Make the same for custom configured plots

# TODO: Make single plot and multi plot decorator that can the quickly execute the plots
# as well as store data in the class about the plot info

# TODO: Make automatic plotting functionality with axarr -> take lengths for components

# The data path to the general data
DATA_PATH = "/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142666/PRODUCT/UTs/20190514/"

# Different baseline-configurations (small-, medium-, large) AT & UT.
# Telescope names and "sta_index"
SMALL = {1: "A0", 5: "B2", 13: "D0", 10: "C1"}
MED = {28: "K0", 18: "G1", 13: "D0", 24: "J3"}
LARGE = {1: "A0", 18: "G1", 23: "J2", 24: "J3"}
UT = {32: "UT1", 33: "UT2", 34: "UT3", 35: "UT4"}

# All baseline configurations
ALL_TELS = {}
ALL_TELS.update(SMALL)
ALL_TELS.update(MED)
ALL_TELS.update(LARGE)
ALL_TELS.update(UT)


class Plotter:
    """Class that plots models as well as vis-, t3phi- and uv-data

    Attributes
    ----------
    """
    def __init__(self, fits_files: List,
                 flux_file: Optional[Path] = "",
                 axis_cutoffs: Optional[List] = [11, -17],
                 limit_spacing: Optional[float] = 0.05,
                 save_path: Optional[Path] = "") -> None:
        """Initialises the class instance"""
        self.fits_files = fits_files

        if len(self.fits_files) > 1:
            warn("There is more than one (.fits)-file detected."\
                 " WARNING: Only data of SAME BAND can be unified into one plot!")

        self.cutoff_low, self.cutoff_high = axis_cutoffs
        self.limit_spacing = limit_spacing

        self.all_vis, self.all_vis_err = [], []
        self.all_vis2, self.all_vis2_err = [], []
        self.all_t3phi, self.all_t3phi_err = [], []
        self.all_flux, self.all_flux_err = [], []

        # TODO: Properly implement different baseline configurations
        self.all_tel_vis, self.all_tel_t3phi = [], []
        self.all_ucoords, self.all_vcoords = [], []

        self.fits_files_names = []

        for fits_file in self.fits_files:
            self.readout = ReadoutFits(fits_file)
            self.fits_file_name = os.path.basename(fits_file)
            self.fits_files_names.append(self.fits_file_name)
            # TODO: Make proper fits-file names for flux plots

            # Fetches all the relevant data from the '.fits'-file
            self.visdata, self.viserr = map(lambda x: x[:6],
                                              self.readout.get_vis()[:2])
            self.vis2data, self.vis2err = map(lambda x: x[:6],
                                              self.readout.get_vis2()[:2])
            self.t3phidata, self.t3phierr = map(lambda x: x[:4],
                                                self.readout.get_t3phi()[:2])

            if not flux_file:
                try:
                    self.flux = self.readout.get_flux()
                except:
                    self.flux = None
                    print("No total flux found!")
            else:
                # FIXME: Make get data function like for fitting
                self.flux = flux_file
            if self.flux is not None:
                self.all_flux.append(self.flux)

            # Get the telescope indices
            self.vissta = self.readout.get_vis()[2]
            self.t3phista = self.readout.get_t3phi()[2]

            # Sets the descriptors of the telescopes' baselines and the closure # phases
            self.tel_vis = np.array([("-".join([ALL_TELS[t]\
                                         for t in duo])) for duo in self.vissta])
            self.tel_t3phi = np.array([("-".join([ALL_TELS[t] for t in trio]))\
                                       for trio in self.t3phista])

            # Gets other important data
            self.ucoords, self.vcoords = map(lambda x: x[:6],\
                                             self.readout.get_split_uvcoords())
            self.wl = self.readout.get_wl()

            self.all_vis.extend(self.visdata)
            self.all_vis_err.extend(self.viserr)
            self.all_vis2.extend(self.vis2data)
            self.all_vis2_err.extend(self.vis2err)
            self.all_t3phi.extend(self.t3phidata)
            self.all_t3phi_err.extend(self.t3phierr)

            self.all_tel_vis.extend(self.tel_vis)
            self.all_tel_t3phi.extend(self.tel_t3phi)

            self.all_ucoords.extend(self.ucoords)
            self.all_vcoords.extend(self.vcoords)

            # The mean of the wavelength.
            # The mean of all the visibilities and their standard deviation
            # self.mean_wl = np.mean(self.wl)
            # self.wl_slice= [j for j in self.wl\
                            # if (j >= self.mean_wl-0.5e-06 and j <= self.mean_wl+0.5e-06)]
            # self.si, self.ei = (int(np.where(self.wl == self.wl_slice[0])[0])-5,
                                    # int(np.where(self.wl == self.wl_slice[~0])[0])+5)

            # self.mean_bin_vis2 = [np.nanmean(i[self.si:self.ei]) for i in self.vis2data]
            # self.baseline_distances = [np.sqrt(x**2+y**2)\
                                       # for x, y in zip(self.ucoords, self.vcoords)]
        self.wl = self.wl[self.cutoff_low:self.cutoff_high]
        self.all_flux = [i[self.cutoff_low:self.cutoff_high] for i in self.all_flux]
        self.all_vis = [i[self.cutoff_low: self.cutoff_high] for i in self.all_vis]
        self.all_vis2 = [i[self.cutoff_low: self.cutoff_high] for i in self.all_vis2]
        self.all_t3phi = [i[self.cutoff_low: self.cutoff_high] for i in self.all_t3phi]

    def get_plot_lims(self, data: List):
        """Gets the low and high limit for a plot from its data

        data: List
            The data from which the limits are to be determined
        """
        low_lim, high_lim = np.min(data), np.max(data)
        y_space = np.sqrt(low_lim**2+high_lim**2)*self.limit_spacing
        return [low_lim-y_space, high_lim+y_space]

    def get_plot_linestyle(self, already_chosen_linestyles: List):
        """Gets a linestyle, which is different from the already chosen one

        Parameters
        ----------
        already_chosen_linestyles: List
            A list that contains the linestyles that have already been applied
        """
        linestyles = ["-", "--", "-.", ":"]
        for linestyle in linestyles:
            if linestyle not in already_chosen_linestyles:
                return linestyle

        if all(linestyle in linestyles in already_chosen_linestyles):
            raise RuntimeError("No linestyle to pick as too many have been picked already!")

    def plot_data(self, ax, data: List, data_name: str,
                  tel_data: Optional[List] = [],
                  lband: Optional[bool] = False,
                  unwrap: Optional[bool] = False) -> None:
        """A template to plot multiple or a singular datasets

        Parameters
        ----------
        ax
            The matplotlib axis to be used for the plot
        data: List
            The plots dataset
        data_name: str
            The name of the data to be plotted. Determines the legend- and plot labels
        tel_data: List, optional
            The telescops' names
        lband: bool, optional
            This specifies the axis used for cutoff, if set to "True" then region between
            L&M band, which has high peaks, is ignored
        unwrap: bool, optional
            This can only be activated if "data_name == cphase". It unwraps the phase from
            180 deg for better visualisation
        """
        already_chosen_linestyles = []
        for i, dataset in enumerate(data):
            if (data_name == "vis") or (data_name == "vis2")\
               or (data_name == "corr_flux"):
                if not tel_data:
                    raise IOError("For this dataset 'tel_data' is needed!")
                baseline = np.around(np.sqrt(self.all_ucoords[i]**2+\
                                             self.all_vcoords[i]**2), 2)
                pas = np.around((np.degrees(
                    np.arctan2(self.all_vcoords[i], self.all_ucoords[i]))-90)*-1, 2)
                label = fr"{tel_data[i]} $B_p$={baseline} m $\phi={pas}^\circ$"
            elif data_name == "cphases":
                if not tel_data:
                    raise IOError("For this dataset 'tel_data' is needed!")
                label = tel_data[i]
            elif data_name == "flux":
                label = None
            else:
                raise IOError(f"The 'data_name': {data_name} does not exist!")

            if data_name == "cphases":
                if unwrap:
                    dataset = np.unwrap(dataset, period=180)
                if i % 4 == 0:
                    linestyle = self.get_plot_linestyle(already_chosen_linestyles)
                    already_chosen_linestyles.append(linestyle)
            elif data_name == "flux":
                linestyle = "-"
                already_chosen_linestyles.append(linestyle)
            else:
                if i % 6 == 0:
                    linestyle = self.get_plot_linestyle(already_chosen_linestyles)
                    already_chosen_linestyles.append(linestyle)

            ax.plot(self.wl*1e6, dataset, linestyle=linestyle, label=label, linewidth=2)

        if data_name == "vis":
            ax.set_ylabel("Vis")
        elif data_name == "vis2":
            ax.set_ylabel("Vis2")
        elif data_name == "corr_flux":
            ax.set_ylabel("Corr. flux [Jy]")
        elif data_name == "flux":
            ax.set_ylabel("Flux [Jy]")
        else:
            if unwrap:
                ax.set_ylabel("Closure phases - Unwrapped [deg]")
            else:
                ax.set_ylabel("Closure phases [deg]")

        if lband:
            indices = (np.where(self.wl < 4.25e-6)[0].tolist() +\
                       np.where(self.wl > 4.45e-6)[0].tolist()).sort()
            plot_lims = self.get_plot_lims([i[0][indices] for i in data.copy()])
        else:
            plot_lims = self.get_plot_lims(data)
        ncol = len(data)//4 if data_name == "cphases" else len(data)//6

        ax.set_xlabel(r"Wavelength [$\mathrm{\mu}$m]")
        ax.set_ylim(plot_lims)
        ax.legend(loc=2, prop={'size': 6}, ncol=ncol)

    def plot_flux(self, ax) -> None:
        """Plots the flux

        Parameters
        ----------
        ax
            The matplotlib axis to be used for the plot
        """
        self.plot_data(ax, self.all_flux, "flux")

    def plot_vis(self, ax, lband: Optional[bool] = False) -> None:
        """Plots all the visibilities/correlated_fluxes

        Parameters
        ----------
        ax
            The matplotlib axis to be used for the plot
        lband: bool, optional
            This specifies the axis used for cutoff, if set to "True" then region between
            L&M band, which has high peaks, is ignored
        """
        self.plot_data(ax, self.all_vis, "vis", self.all_tel_vis, lband)

    def plot_corr_flux(self, ax, lband: Optional[bool] = False) -> None:
        """Plots all the visibilities/correlated_fluxes in one plot

        Parameters
        ----------
        ax
            The matplotlib axis to be used for the plot
        lband: bool, optional
            This specifies the axis used for cutoff, if set to "True" then region between
            L&M band, which has high peaks, is ignored
        """
        self.plot_data(ax, self.all_vis, "corr_flux", self.all_tel_vis, lband)

    def plot_vis2(self, ax, lband: Optional[bool] = False) -> None:
        """Plots all the visibilities squared in one plot

        Parameters
        ----------
        ax
            The matplotlib axis to be used for the plot
        """
        self.plot_data(ax, self.all_vis2, "vis2", self.all_tel_vis, lband)

    def plot_cphases(self, ax,
                     lband: Optional[bool] = False,
                     unwrap: Optional[bool] = False) -> None:
        """Plots all the closure phases into one plot

        Parameters
        ----------
        ax
            The matplotlib axis to be used for the plot
        lband: bool, optional
            This specifies the axis used for cutoff, if set to "True" then region between
            L&M band, which has high peaks, is ignored
        unwrap: bool, optional
            This can only be activated if "data_name == cphase". It unwraps the phase from
            180 deg for better visualisation
        """
        self.plot_data(ax, self.all_t3phi, "cphases", self.all_tel_t3phi, lband, unwrap)

    def plot_vis24baseline(self, ax, do_fit: Optional[bool] = True) -> None:
        """ Plot the mean visibility for one certain wavelength and fit it with a
        gaussian and airy disk

        Parameters
        ----------
        ax
            The matplotlib axis to be used for the plot
        do_fit: bool, optional
            Switches the fit mechanic on or off. DEFAULT "True"
        """
        ax.plot(self.baseline_distances, self.mean_bin_vis2, ls='None', marker='o')

        # TODO: Implement fit here
        if do_fit:
            ...

        ax.set_xlabel(fr'uv-distance [m] at $\lambda_0$={self.mean_wl:.2f} $\mu m$')
        ax.set_ylabel(r'$\bar{V}$')

    def plot_waterfall(self, ax) -> None:
        """Plots a waterfall with the mean wavelength for the different baselines"""
        for i in range(6):
            ax.errorbar(self.wl[self.si:self.ei]*1e06,
                        self.vis2data[i][self.si:self.ei],
                        yerr=np.nanstd(self.vis2data[i][self.si:self.ei]),
                        label=self.tel_vis2[i], ls='None', fmt='o')
            ax.set_xlabel(r'wl [micron]')
            ax.set_ylabel('vis2')
            ax.legend(loc='best')

    def plot_uv(self, ax) -> None:
        """Plots the uv-coordinates with an orientational compass and the synthethis

        Parameters
        ----------
        ax
            The axis anchor of matplotlib.pyplot

        Returns
        -------
        None
        """
        ax.scatter(self.ucoords, self.vcoords)
        ax.scatter(-self.ucoords, -self.vcoords)
        ax.set_xlim([150, -150])
        ax.set_ylim([-150, 150])
        ax.set_ylabel('v [m]')
        ax.set_xlabel('u [m]')

        # Compass for the directions
        cardinal_vectors = [(0,1), (0,-1), (1,0), (-1,0)]
        cardinal_colors  = ['black', 'green', 'blue', 'red']
        cardinal_directions = ['N', 'S', 'W', 'E']
        arrow_size, head_size = 40, 10
        x, y = (-85, 85)

        for vector, color, direction in zip(cardinal_vectors,
                                            cardinal_colors, cardinal_directions):
            dx, dy = vector[0]*arrow_size, vector[1]*arrow_size
            if vector[0] == 0:
                ax.text(x-dx-5, y+dy, direction)
            if vector[1] == 0:
                ax.text(x-dx, y+dy+5, direction)
            arrow_args = {"length_includes_head": True,
                          "head_width": head_size, "head_length": head_size,
                          "width": 1, "fc": color, "ec": color}
            ax.arrow(x, y, dx, dy, **arrow_args)


def rotation_synthesis_uv(inp):
    """This function was written by Jozsef Varga (from menEWS: menEWS_plot.py).

    Calculates uv-point corresponding to inp (see "get_header_info"),
    for hour angle(s) (ha)
    """
    ra, dec, BE, BN, BL, base = inp
    paranal_lat = -24.62587 * np.pi / 180.

    u = BE * np.cos(ha) -\
            BN * np.sin(lat) * np.sin(ha) + BL * np.cos(lat) * np.sin(ha)
    v = BE * np.sin(dec) * np.sin(ha) +\
            BN * (np.sin(lat) * np.sin(dec) * np.cos(ha) +\
                  np.cos(lat) * np.cos(dec)) - BL * \
        (np.cos(lat) * np.sin(dec) * np.cos(ha)- np.sin(lat) * np.cos(dec))
    return u, v


def make_uv_tracks(uv, inp, flag, ax, bases=[], symbol='x',color='',
    print_station_names=True,sel_wl=1.0,plot_Mlambda=False):
    """This function was written by Jozsef Varga (from menEWS: menEWS_plot.py).

    From coordinate + ha (range), calculate uv tracks"""

    ra, dec, BE, BN, BL, base = inp
    paranal_lat = -24.62587 * np.pi / 180.
    mlim = 2.0  # airmass limit for tracks

    if plot_Mlambda == True:
        u, v = map(lambda x: x/sel_wl, uv)
    else:
        u, v = uv

    if not color:
        if np.all(flag) == 'True':
            color = 'r'
        else:
            color = 'g'

    if base not in bases:
        hamax = np.arccos(abs((1. / mlim - np.sin(lat) * np.sin(dec)) / \
                              (np.cos(lat) * np.cos(dec))))
        harng = np.linspace(-hamax, hamax, 1000)

        ul, vl = ulvl = calculate_uv_points(inp, harng)
        if plot_Mlambda == True:
            u, v = map(lambda x: x/sel_wl, ulvl)

        ax.plot(ul, vl, '-', color='grey',alpha=0.5)
        ax.plot(-ul, -vl, '-', color='grey',alpha=0.5)
        ax.plot([0.], [0.], '+k', markersize=5, markeredgewidth=2,alpha=0.5)

        if print_station_names:
            ax.text(-u-7, -v-3, base, color='0',alpha=0.8)
        bases.append(base)

    ax.plot(u, v, symbol, color=color, markersize=10, markeredgewidth=3)
    ax.plot(-u, -v, symbol, color=color, markersize=10, markeredgewidth=3)

    return bases


def make_uv_plot(dic,ax,verbose=False,annotate=True,B_lim=(np.nan,np.nan),figsize=(5,5),
    color='',print_station_names=True,sel_wl=1.0,plot_Mlambda=False):
    """This function was written by Jozsef Varga (from menEWS: menEWS_plot.py)"""
    if plot_Mlambda==False:
        sel_wl = 1.0
    try:
        u = dic['VIS2']['U']
        v = dic['VIS2']['V']
        flag = dic['VIS2']['FLAG']
        sta_index = dic['VIS2']['STA_INDEX']
        mjd = dic['VIS2']['MJD']
    except KeyError as e:
        if verbose: print(e)
        u = [0.0]
        v = [0.0]
        flags = [False]
        sta_index = []
        mjd = [0.0]

    uvs = []
    inps = []
    flags = []
    umax = []
    vmax = []
    for j in range(len(u)):
        uvs.append([u[j],v[j]])
        try:
            BE, BN, BL = dic['STAXYZ'][sta_index[j, 0] == dic['STA_INDEX']][0] - \
                dic['STAXYZ'][sta_index[j, 1] == dic['STA_INDEX']][0]
            sta_label= dic['STA_NAME'][sta_index[j, 0] == dic['STA_INDEX']][0] + '-' + \
                        dic['STA_NAME'][sta_index[j, 1] == dic['STA_INDEX']][0]
        except IndexError as e:
            print('make_uv_plot STA_INDEX error.')
            print(e)
            BE, BN, BL = [np.nan,np.nan,np.nan]
            sta_label= ''
        inps.append( [dic['RA'] * np.pi / 180., dic['DEC'] * np.pi / 180., BE, BN, BL, sta_label]  )
        flags.append(flag[j])
    bases = []
    umax = np.nanmax(np.abs(u))
    vmax = np.nanmax(np.abs(v))
    if not (dic['MJD-OBS']):
        dic['MJD-OBS'] = np.amin(mjd[0])
    try:
        rel_time = (mjd - dic['MJD-OBS']) * 24.0 * 3600.0  # (s)
        dic['TREL'] = rel_time[0]

        for k, uv in enumerate(uvs):
            bases = make_uv_tracks(uv, inps[k], flags[k],ax, bases,
            color=color,print_station_names=print_station_names,
            sel_wl=sel_wl,plot_Mlambda=plot_Mlambda)

        if plot_Mlambda == False:
            xlabel ='$u$ (m)'
            ylabel ='$v$ (m)'
        else:
            xlabel ='$u$ ($M\lambda$)'
            ylabel ='$v$ ($M\lambda$)'
        ax.set_xlim((130, -130))
        ax.set_ylim((-130, 130))
        plotmax = 1.3*np.amax([umax,vmax])

        plot_title = dic['TARGET'] + "\n" + "date: " + dic['DATE-OBS'] + "\n" + "TPL start: " + dic['TPL_START'] + "\n" + dic['CATEGORY'] + ' ' +\
            dic['BAND'] + ' ' + dic['DISPNAME'] #+ ' ' + dic['BCD1'] + '-' + dic['BCD2']
        if math.isnan(B_lim[0]):
            xlim = (+plotmax/ sel_wl,-plotmax/ sel_wl)
            ylim = (-plotmax/ sel_wl,+plotmax/ sel_wl)
        else:
            xlim = (+B_lim[1]/ sel_wl,-B_lim[1]/ sel_wl)
            ylim = (-B_lim[1]/ sel_wl,+B_lim[1]/ sel_wl)
        #if plot_Mlambda == True:
        plot_config(xlabel, ylabel,plot_title, ax, dic,
                    ylim=ylim,xlim=xlim,plot_legend=False,annotate=annotate)
    except TypeError as e:
        if verbose: print('Unable to plot ' + 'uv')
        if verbose: print(e)
        return 1

    return 0


if __name__ == ('__main__'):
    specific_path = "combined_and_averaged/split"
    fits_files_N = ["2019-05-14T05_28_03.2019-05-14T04_52_11_cal_avg_N_cphases.fits"]
    fits_files_N = [os.path.join(DATA_PATH, specific_path, i) for i in fits_files_N]
    plotter_N = Plotter(fits_files_N)
    fig, axarr = plt.subplots(1, 2)
    ax, bx = axarr.flatten()
    plotter_N.plot_corr_flux(ax)
    plotter_N.plot_cphases(bx)
    plt.savefig("test.png")

