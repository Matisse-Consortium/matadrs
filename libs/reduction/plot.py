import os
from pathlib import Path
from warnings import warn
from typing import List, Optional

import numpy as np
import matplotlib.pyplot as plt

from matadrs.libs.reduction.readout import ReadoutFits

# TODO: Make either model or fourier transform carry more info like the name of
# the plot or similar -> Work more with classes

# TODO: Implement plot save func

# TODO: Make standalone method that makes plot plot immediately
# TODO: Make the same for custom configured plots

# TODO: Make single plot and multi plot decorator that can the quickly execute the plots
# as well as store data in the class about the plot info

# TODO: Make automatic plotting functionality with axarr -> take lengths for components

# The data path to the general data
DATA_DIR = Path(__file__).parent.parent.parent.parent.parent / "data"


class Plotter:
    """Class that plots models as well as vis-, t3phi- and uv-data

    Attributes
    ----------
    """
    def __init__(self, fits_files: List[Path],
                 flux_file: Optional[Path] = None,
                 lband: Optional[bool]= False,
                 limit_spacing: Optional[float] = 0.05,
                 plot_name: Optional[str] = None,
                 axis_cutoffs: Optional[List] = [0, 0],
                 save_path: Optional[Path] = None) -> None:
        """Initialises the class instance"""
        self.lband = lband

        if plot_name is None:
            self.plot_name = "combined_files.pdf" if len(fits_files) > 1\
                else f"{os.path.basename(fits_files[0]).split('.')[0]}.pdf"
        else:
            self.plot_name = plot_name

        self.save_path = save_path
        self.components = []

        if len(fits_files) > 1:
            warn("There is more than one (.fits)-file detected."\
                 " WARNING: Only data of SAME BAND can be unified into one plot!")

        self.cutoff_low, self.cutoff_high = axis_cutoffs
        self.limit_low = self.cutoff_low if self.cutoff_low > 13 else 13
        self.limit_high = self.cutoff_high if self.cutoff_high < -17 else -17
        self.limit_spacing = limit_spacing

        self.readout = ReadoutFits(fits_files)

    @property
    def number_of_plots(self):
        return len(self.components)

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
            raise RuntimeError("No linestyle to pick as all have been picked already!")

    def plot_data(self, ax, data: List, data_name: str,
                  tel_data: Optional[List] = [],
                  lband: Optional[bool] = False,) -> None:
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
        """
        already_chosen_linestyles = []
        for index, dataset in enumerate(data):
            if (data_name == "vis") or (data_name == "vis2")\
               or (data_name == "corr_flux"):
                # TODO: Make this with new readout
                baseline = np.around(np.sqrt(self.all_ucoords[index]**2+\
                                             self.all_vcoords[index]**2), 2)
                pas = np.around((np.degrees(np.arctan2(self.readout.oi_vis[index],
                                                       self.all_ucoords[index]))-90)*-1, 2)
                label = fr"{tel_data[index]} $B_p$={baseline} m $\phi={pas}^\circ$"
            else:
                label = tel_data[index]

            if data_name == "cphases":
                if index % 4 == 0:
                    linestyle = self.get_plot_linestyle(already_chosen_linestyles)
                    already_chosen_linestyles.append(linestyle)
            elif data_name == "flux":
                linestyle = "-"
                already_chosen_linestyles.append(linestyle)
            else:
                if index % 6 == 0:
                    linestyle = self.get_plot_linestyle(already_chosen_linestyles)
                    already_chosen_linestyles.append(linestyle)

            ax.plot(self.wl*1e6, dataset, linestyle=linestyle, label=label, linewidth=2)

        if data_name == "vis":
            ax.set_ylabel("Vis")
        elif data_name == "vis2":
            ax.set_ylabel("Vis2")
        elif data_name == "corr_flux":
            ax.set_ylabel("Correlated flux [Jy]")
        elif data_name == "flux":
            ax.set_ylabel("Flux [Jy]")
        else:
            ax.set_ylabel(r"Closure phases [$^{\circ}$]")

        if data_name == "flux":
            data_bounding = [dataset[:self.limit_high] for dataset in data]
        else:
            data_bounding = [dataset[self.limit_low:self.limit_high]\
                             for dataset in data]

        # TODO: Fix this functionality -> Set the limits properly
        # TODO: Set no limits for the closure phases
        if lband:
            wavelength = self.wl[self.limit_low:self.limit_high]*1e6
            lower = np.where(wavelength < 4.2)[0].tolist()
            upper = np.where(wavelength > 4.6)[0].tolist()
            indices = upper + lower
            plot_lims = self.get_plot_lims([dataset[indices]\
                                            for dataset in data_bounding])
        else:
            plot_lims = self.get_plot_lims(data_bounding)

        plot_lims_after = [limit+0.25*limit for limit in plot_lims]
        ncol = len(data)//4 if data_name == "cphases" else\
                (len(data)//6 if data_name in\
                 ["vis", "vis2", "corr_flux"] else 1)

        ax.set_xlabel(r"Wavelength [$\mathrm{\mu}$m]")
        ax.set_ylim(plot_lims)
        ax.legend(loc=1, prop={'size': 6}, ncol=ncol)

    def _make_component(self, data: List, data_name: str,
                        tel_data: Optional[List] = [],
                        lband: Optional[bool] = False):
        """Makes a plot component

        Parameters
        ----------
        data: List
            The plots dataset
        data_name: str
            The name of the data to be plotted. Determines the legend- and plot labels
        tel_data: List, optional
            The telescops' names
        lband: bool, optional
            This specifies the axis used for cutoff, if set to "True" then region between
            L&M band, which has high peaks, is ignored

        Returns
        -------
        Dict
            The component's information
        """
        return {"data": data, "data_name": data_name,
                "tel_data": tel_data, "lband": lband}

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
            for ax, component in zip(axarr.flatten(), self.components):
                self.plot_data(ax, *component.values())
        else:
            self.plot_data(axarr, *self.components[0].values())

        fig.tight_layout()

        if save:
            plt.savefig(os.path.join(self.save_path, self.plot_name), format="pdf")
        else:
            plt.show()
        plt.close()

    def add_flux(self) -> None:
        """Plots the flux """
        self.components.append(self._make_component(self.all_flux, "flux",
                                                    self.all_tel_flux))
        return self

    def add_vis(self) -> None:
        """Plots all the visibilities/correlated_fluxes """
        self.components.append(self._make_component(self.all_vis, "vis",
                                                    self.all_tel_vis,
                                                    self.lband))
        return self

    def add_corr_flux(self, lband: Optional[bool] = False) -> None:
        """Plots all the visibilities/correlated_fluxes in one plot """
        self.components.append(self._make_component(self.all_vis, "corr_flux",
                                                    self.all_tel_vis,
                                                    self.lband))
        return self

    def add_vis2(self) -> None:
        """Plots all the visibilities squared in one plot"""
        self.components.append(_make_component(self.all_vis2, "vis2",
                                               self.all_tel_vis,
                                               self.lband))
        return self

    def add_cphases(self) -> None:
        """Plots all the closure phases into one plot"""
        self.components.append(self._make_component(self.all_t3phi, "cphases",
                                                    self.all_tel_t3phi,
                                                    self.lband))
        return self

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
    # TODO: Fix the plotter for the pandas dataframe structure
    fits_files = DATA_DIR / "jozsef_reductions" / "hd142666_2019-05-14T05_28_03_N_TARGET_FINAL_INT.fits"
    plot_fits = Plotter([fits_files])
    plot_fits.add_cphases().add_corr_flux().add_flux().plot()


