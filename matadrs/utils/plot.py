import re
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Callable, Tuple, List, Union, Optional

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from .readout import ReadoutFits
from .tools import unwrap_phases, calculate_uv_tracks, get_fits_by_tag, \
        get_colorlist
from .options import OPTIONS
from ..mat_tools.mat_show_atmo_param_v2 import show_seeing
from ..mat_tools.mat_show_oifits_pbe_ama_short import open_oi_dir, \
        filter_oi_list, show_vis_tf_vs_time


matplotlib.use('Agg')


def plot_data_quality(
        reduced_directory: Path, output_dir: Path) -> None:
    """Plots the data quality of the reduced fits-files."""
    if not output_dir.exists():
        output_dir.mkdir(parents=True)

    fits_file = get_fits_by_tag(reduced_directory, "RAW_INT")
    if not fits_file:
        return
    readout = ReadoutFits(fits_file[0])

    dics = open_oi_dir(reduced_directory, choice_band_LM="L")
    res = readout.resolution.upper()
    plot_kwargs = {"output_path": output_dir, "date": readout.date,
                   "target": readout.name, "wlenRange": [3.2, 3.9],
                   "saveplots": True, "show": False, "plot_errorbars": False}

    dics = filter_oi_list(
            dics, spectral_resolutions=[res],
            DIT_range=[0.111, 11.], dates=[readout.date], bands=["L"])

    show_seeing(dics, **plot_kwargs)
    show_vis_tf_vs_time(dics, **plot_kwargs)


def plot_broken_axis(ax: Axes, x: np.ndarray,
                     y: np.ndarray, yerr: np.ndarray,
                     range1: Tuple[float, float],
                     range2: Tuple[float, float],
                     ax_left: Optional[Axes] = None,
                     ax_right: Optional[Axes] = None,
                     color: Optional[str] = None,
                     ylims: Optional[Tuple[float, float]] = None,
                     error: Optional[bool] = False,
                     no_fill: Optional[bool] = False,
                     err_percentile: Optional[float] = None,
                     **kwargs):
    """Plot two axes next to each other to display the L-band.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The original axis to be split.
    x : np.ndarray
        The x-coordinates of the data.
    y : np.ndarray
        The y-coordinates of the data.
    yerr : np.ndarray
        The y-error of the data.
    range1 : tuple
        The range of the first axis.
    range2 : tuple
        The range of the second axis.
    ax_right : matplotlib.axes.Axes, optional
        The left axis.
    ax_right : matplotlib.axes.Axes, optional
        The right axis.
    color : str, optional
        The color of the line.
    ylims : tuple, optional
        The limits of the y-axis.
    label : str, optional
        The label of the first axis.
    error : bool, optional
        If the error should be plotted.
    no_fill : bool, optional
        If the fill between the errors should not be plotted and instead
        be indicated by one errorbar.
    err_percentile : int, optional
        The percentage at which point the sole errorbar will be plotted.
    """
    x1 = x[(x >= range1[0]) & (x <= range1[1])]
    y1 = y[(x >= range1[0]) & (x <= range1[1])]
    yerr1 = yerr[(x >= range1[0]) & (x <= range1[1])]

    x2 = x[(x >= range2[0]) & (x <= range2[1])]
    y2 = y[(x >= range2[0]) & (x <= range2[1])]
    yerr2 = yerr[(x >= range2[0]) & (x <= range2[1])]

    if ax_left is None and ax_right is None:
        ax_left = inset_axes(
                ax, width="48%", height="100%", loc='center left')
        ax_right = inset_axes(
                ax, width="48%", height="100%", loc='center right')

        ax_left.spines.right.set_visible(False)
        ax_right.spines.left.set_visible(False)

        ax_left.yaxis.tick_left()
        ax_left.tick_params(labelright='off')
        ax_right.yaxis.set_visible(False)

        left_ticks = np.around(np.linspace(range1[0], range1[1], 5)[:-1], 1)
        right_ticks = np.around(np.linspace(range2[0], range2[1], 5)[1:], 1)

        ax_left.set_xticks(left_ticks)
        ax_right.set_xticks(right_ticks)

        ax_left.set_ylim(*ylims)
        ax_right.set_ylim(*ylims)

        ax_left.set_xlim(range1[0], range1[1])
        ax_right.set_xlim(range2[0], range2[1])

    ax_left.plot(x1, y1, color=color, **kwargs)
    ax_right.plot(x2, y2, color=color, **kwargs)
    if error:
        if no_fill:
            err_percentile = 0.25 if err_percentile is None else err_percentile
            no_fill_index = int(np.ceil(x1.shape[-1]*err_percentile))
            ax_left.errorbar(
                x1[no_fill_index], y1[no_fill_index],
                yerr=np.mean(yerr1), color=color, capsize=3)
        else:
            ax_left.fill_between(x1, y1+yerr1, y1-yerr1,
                                 color=color, alpha=0.2)
            ax_right.fill_between(x2, y2+yerr2, y2-yerr2,
                                  color=color, alpha=0.2)

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    return ax_left, ax_right


@dataclass
class PlotComponent:
    """Class containing the elements required for a plot
    from the reduced data.

    Parameters
    ----------
    labels : list
    x_values : list
    y_values : list
    y_errors : list
    band : str
    """
    labels: List = None
    x_values: List = None
    y_values: List = None
    y_errors: List = None


# TODO: Implement text plotter with the information on the observation
class Plotter:
    """Class that plots models as well as reduced data

    Parameters
    ----------
    fits_files : pathlib.Path or list of pathlib.Path
    plot_name : str, optional
    save_path : pathlib.Path, optional

    Attributes
    ----------
    num_components : int
        The number of components in the plot.

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
    add_vis(legend_format="long")
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

    def __init__(self, fits_files: Union[Path, List[Path]],
                 plot_name: Optional[str] = None,
                 save_dir: Optional[Path] = None) -> None:
        """The class's constructor"""
        self.fits_files = [fits_files]\
            if not isinstance(fits_files, List) else fits_files

        if save_dir is None:
            self.save_dir = Path("").cwd()
        else:
            self.save_dir = Path(save_dir)

        if plot_name is not None:
            self.plot_name = plot_name
        elif isinstance(fits_files, List) and len(fits_files) > 1:
            self.plot_name = "combined_fits.pdf"
        else:
            self.plot_name = f"{Path(fits_files).stem}.pdf"

        self.color_grouping = "file"
        self.lband_bounds = (3.1, 3.9)
        self.mband_bounds = (4.55, 4.9)
        self.nband_bounds = (8.5, 12.5)

        self.readouts_unfiltered = list(map(ReadoutFits, self.fits_files))
        self.readouts = self.readouts_unfiltered.copy()
        self.sort("date")
        self.components = {}

    def __str__(self):
        """The class's string representation"""
        return f"Plotting the following (.fits)-files:\n{'':-^50}\n" + \
            "\n".join([readout.fits_file.stem for readout in self.readouts])

    @property
    def num_components(self):
        """The number of componets contained"""
        return len(self.components)

    def _set_y_limits(self, wavelength: np.ndarray,
                      data: List[np.ndarray],
                      margin: Optional[float] = 0.05) -> Tuple[int, int]:
        """Sets the y-limits from the data with some margin"""
        try:
            if np.min(wavelength) >= 6:
                indices = np.where((wavelength > self.nband_bounds[0])
                                   | (wavelength < self.nband_bounds[1]))
            else:
                indices_high = np.where((wavelength >= self.mband_bounds[0])
                                       & (wavelength <= self.mband_bounds[1]))
                indices_low = np.where((wavelength >= self.lband_bounds[0])
                                        & (wavelength <= self.lband_bounds[1]))
                indices = np.hstack((indices_high, indices_low))
            ymin, ymax = data[:, indices].min(), data[:, indices].max()
        except ValueError:
            ymin, ymax = np.percentile(data, 10), np.percentile(data, 90)
        spacing = np.linalg.norm(ymax-ymin)*margin
        if np.isnan(ymin) and np.isnan(ymax):
            return None, None
        return ymin-spacing, ymax+spacing

    def sort(self, by: str):
        """Filters the object's data by the given key.

        The key can be any of the properties of the ReadoutFits class

        Parameters
        ----------
        by : str
            The data to sort by.
        """
        self.readouts = sorted(self.readouts, key=lambda x: getattr(x, by.lower()))
        return self

    def filter(self, by: List[str], contains: List[str]):
        """Filters the object's data by the given key.

        The key can be any of the properties of the ReadoutFits class

        Parameters
        ----------
        by : str
            The key to filter the data by.
        contains : str
            The string that the key should contain.
        """
        by = [by] if not isinstance(by, List) else by
        contains = [contains] if not isinstance(contains, List) else contains

        for key, contain in zip(by, contains):
            self.readouts = [readout for readout in self.readouts if contain.lower() \
                in getattr(readout, key.lower()).lower()]
        return self

    def plot_uv(self, ax: Axes,
                symbol: Optional[str] = "x",
                airmass_lim: Optional[float] = 2.,
                show_text: Optional[List] = False,
                make_tracks: Optional[bool] = True,
                show_legend: Optional[bool] = True,
                legend_location: Optional[str] = OPTIONS.plot.legend.location,
                legend_size: Optional[int] = OPTIONS.plot.legend.fontsize,
                color_by: Optional[str] = "file",
                uv_extent: Optional[int] = None,
                readouts: Optional[ReadoutFits] = None,
                **kwargs) -> None:
        """Plots the (u, v)-coordinates and their corresponding tracks

        Parameters
        ----------
        ax : matplotlib.axes.Axes
        symbol : str, optional
            The symbol that markes the coordinates of the (u, v)-point.
        airmass_lim: float, optional
            The airmass limit for the uv-coords.
        show_text: list of bool, optional
            If the baselines should be shown next to the coordinates as text.
        make_tracks: bool, optional
            If the tracks should be plotted.
        show_legend: bool, optional
            If the legend should be shown.
        legend_location : str, optional
            The location of the legend.
        legend_size : int, optional
            The size of the legend.
        color_by : str, optional
            The color grouping used for the uv-coords. If 'file' the
            colors are based on the different (.fits)-files. If 'instrument'
            the colors are based on the different instruments (see
            ReadoutFits.instrument).
        """
        instruments, handles, uv_max = [], [], 0
        readouts = self.readouts if readouts is None else readouts
        colors = get_colorlist(OPTIONS.plot.color.colormap, len(readouts) * 6 + 1)
        for index, readout in enumerate(readouts):
            uv_coords = readout.oi_vis2["UVCOORD"]
            if uv_max < (tmp_uv_max := uv_coords.max()):
                uv_max = tmp_uv_max
            sta_indices = readout.oi_vis2["STA_INDEX"]
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

            if color_by == "file":
                color = colors[index]
                handles.append(mlines.Line2D(
                    [], [], color=color, marker="X",
                    linestyle="None", label=readout.date[:-8]))
            elif color_by == "instrument":
                if readout.instrument not in instruments:
                    instruments.append(readout.instrument)
                color = colors[instruments.index(readout.instrument)]

            for uv_index, (u_coords, v_coords) in enumerate(uv_coords):
                ax.plot(u_coords, v_coords, symbol, color=color,
                        markersize=10, markeredgewidth=3)
                ax.plot(-u_coords, -v_coords, symbol,
                        color=color, markersize=10, markeredgewidth=3)

                if show_text:
                    ax.text(-u_coords-3.5, -v_coords-1.5, sta_labels[uv_index],
                            fontsize="small", color='0', alpha=0.8)

                if make_tracks:
                    u_coord_tracks, v_coord_tracks = calculate_uv_tracks(
                        baselines[uv_index], readout.dec*np.pi/180, airmass_lim)
                    ax.plot(u_coord_tracks, v_coord_tracks, '-', color='grey', alpha=0.5)
                    ax.plot(-u_coord_tracks, -v_coord_tracks, '-', color='grey', alpha=0.5)

        ax.plot([0.], [0.], '+k', markersize=5, markeredgewidth=2, alpha=0.5)

        # TODO: Implement check or calculation for the orientations
        xlabel, ylabel = "$u$ (m) - South", "$v$ (m) - East"
        uv_extent = int(uv_max + uv_max*0.25) if uv_extent is None else uv_extent

        if color_by == "instrument":
            handles = []
            for index, instrument in enumerate(instruments):
                color = colors[index]
                handles.append(mlines.Line2D(
                    [], [], color=color, marker="X",
                    linestyle="None", label=instrument.upper()))

        if show_legend:
            ax.legend(handles=handles,
                      loc=legend_location, fontsize=legend_size)

        ax.set_xlim([uv_extent, -uv_extent])
        ax.set_ylim([-uv_extent, uv_extent])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

    def make_component(self, data_name: str,
                       legend_format: Optional[str] = "verbose",
                       unwrap: Optional[bool] = False,
                       period: Optional[int] = 360, **kwargs):
        """Generates a pandas DataFrame that has all the plots' information

        Parameters
        ----------
        data_name : str
            The name of the data to be plotted. Determines the legend- and plot
            labels.
        legend_format : str, optional
            Sets the format of the legend: For all information set "verbose" and
            for only station names set "short". Default is "verbose".
        unwrap : bool, optional
            If the phase should be unwrapped.
        period : int, optional
            The unwrapping period.
        """
        component, component_label, component_kwargs = [], None, {**kwargs}
        for readout in self.readouts:
            sub_component = PlotComponent(
                x_values=readout.oi_wavelength["EFF_WAVE"].data.squeeze())
            if data_name == "flux":
                if "FLUXDATA" in readout.oi_flux.columns:
                    sub_component.y_values = readout.oi_flux["FLUXDATA"]
                    component_label = f"Flux ({readout.get_unit('oi_flux', 'fluxdata')})"
                else:
                    sub_component.y_values = readout.oi_flux["FLUX"]
                    component_label = f"Flux ({readout.get_unit('oi_flux', 'flux')})"

                # HACK: To get flux in Jansky, remove if solved
                component_label = f"Flux (Jy)"
                sub_component.y_errors = readout.oi_flux["FLUXERR"]

                if legend_format == "verbose":
                    labels = [f"{readout.instrument.upper()}/{readout.date.split('T')[0]}"]
                else:
                    labels = readout.get_telescopes() if len(sub_component.y_values) > 1 else ["Averaged"]
                sub_component.labels = labels

            elif data_name in ["vis", "vis2", "diff", "corrflux"]:
                station_names = readout.oi_vis2["DELAY_LINE"]
                if legend_format == "verbose":
                    baselines = np.around(readout.oi_vis2["BASELINE"], 2)
                    u_coords = readout.oi_vis2["UVCOORD"][:, 0]
                    v_coords = readout.oi_vis2["UVCOORD"][:, 1]
                    pas = np.around(
                        (np.degrees(np.arctan2(v_coords, u_coords))-90)*-1, 2)

                    # TODO: Make the variables into mathrm
                    labels = [fr"{station_name} $B_p$={baseline} m $\phi={pa}^\circ$"
                              for station_name, baseline, pa in zip(station_names, baselines, pas)]
                else:
                    labels = station_names
                sub_component.labels = labels

                if data_name == "vis":
                    unit = readout.get_unit('oi_vis', 'visamp')
                    sub_component.y_values = readout.oi_vis["VISAMP"]
                    sub_component.y_errors = readout.oi_vis["VISAMPERR"]
                    if unit == "Jy":
                        component_label = "Correlated Fluxes (Jy)"
                    else:
                        component_label = "Visibilities (a.u.)"
                elif data_name == "diff":
                    diff_phases = readout.oi_vis["VISPHI"]
                    diff_phases_err = readout.oi_vis["VISPHIERR"]
                    if unwrap:
                        diff_phases, diff_phases_err = unwrap_phases(
                            diff_phases, diff_phases_err, period)
                    sub_component.y_values = diff_phases
                    sub_component.y_errors = diff_phases_err
                    component_label = r"Differential Phases ($^{\circ}$)"
                elif data_name == "vis2":
                    sub_component.y_values = readout.oi_vis2["VIS2DATA"]
                    sub_component.y_errors = readout.oi_vis2["VIS2ERR"]
                    component_label = "Squared visibilities (a.u.)"
                else:
                    raise KeyError("No data-type of that data name exists!")

            elif data_name == "cphases":
                cphases = readout.oi_t3["T3PHI"]
                cphases_err = readout.oi_t3["T3PHIERR"]
                if unwrap:
                    cphases, cphases_err = unwrap_phases(
                        cphases, cphases_err, period)
                sub_component.labels = readout.oi_t3["TRIANGLE"]
                sub_component.y_values = cphases
                sub_component.y_errors = cphases_err
                component_label = r"Closure Phases ($^{\circ}$)"
            elif data_name == "uv":
                sub_component = self.plot_uv
                component_label = "$(u, v)$-coordinates"
            else:
                raise KeyError("Input data name cannot be queried!")
            component.append(sub_component)
        self.components[component_label] = {"values": component,
                                            "kwargs": component_kwargs}
        return self

    def add_flux(self, **kwargs):
        """Adds the total flux(es) as a subplot."""
        return self.make_component("flux", **kwargs)

    def add_vis(self, **kwargs):
        """Adds the visibilities/correlated fluxes as a subplot."""
        return self.make_component("vis", **kwargs)

    def add_vis2(self, **kwargs):
        """Adds the squared visibilities as a subplot."""
        return self.make_component("vis2", **kwargs)

    def add_t3(self, **kwargs):
        """Adds the closure phases as a subplot."""
        return self.make_component("cphases", **kwargs)

    def add_diff_phases(self, **kwargs):
        """Adds the differential phases as a subplot."""
        return self.make_component("diff", **kwargs)

    def add_uv(self, **kwargs):
        """Adds the (u, v)-coordinates as a subplot."""
        return self.make_component("uv", **kwargs)

    def add_mosaic(self, **kwargs):
        """Combines multiple subplots to produce a mosaic plot."""
        self.add_uv(**{k: v for k, v in kwargs.items() if k != "legend_size"})
        self.add_vis(**kwargs)
        self.add_vis2(**kwargs)
        self.add_flux(**kwargs)
        self.add_diff_phases(**kwargs)
        self.add_t3(**kwargs)
        return self

    # TODO: Add different linestyles for different files and also different colorschemes
    # or rather the option to choose
    def plot_components(
        self, ax: plt.Axes, name: str,
        component: Union[Callable, PlotComponent],
        sharex: Optional[bool] = False,
        show_legend: Optional[bool] = True,
        share_legend: Optional[bool] = False,
        legend_location: Optional[str] = OPTIONS.plot.legend.location,
        legend_size: Optional[int] = OPTIONS.plot.legend.fontsize,
        error: Optional[bool] = False,
        no_fill: Optional[bool] = False,
        margin: Optional[float] = 0.05,
        color_by: Optional[str] = "file",
        **kwargs) -> Axes:
        """Plots all the data of a single component.

        Parameters
        ----------
        fig : matplotlib.figure
            The figure to add the subplot to.
        column : int
            The column of the subplot grid.
        row : int
            The row of the subplot grid.
        name : str
            The name of the component.
        component : callable or plotcomponent
            The component to plot.
        sharex : bool, optional
            If True, the x-axis will be shared.
        show_legend : bool, optional
            If True, the legend will be shown.
        share_legend : bool, optional
            If True, the legend will be shared.
        legend_location : str, optional
            The location of the legend.
        legend_size : int, optional
            The size of the legend.
        error : bool, optional
            If True, the error bars will be plotted.
        no_fill : bool, optional
            If True, the fill between the errors will not be plotted and
            instead be indicated by one errorbar.
        margin : bool, optional
            The margin around the plot.
        color_by : str, optional
        kwargs : dict
        """
        xlabel = r"$\lambda$ ($\mathrm{\mu}$m)"
        colors = get_colorlist(OPTIONS.plot.color.colormap, len(self.readouts) * 6 ** 2)

        # TODO: Make it so that the offsets don't overshoot the lenght of the data
        # TODO: Check if this works still
        offset = 0.8/len(component)
        ax_left, ax_right, handles = None, None, []
        for comp_index, sub_component in enumerate(component, start=1):
            if isinstance(sub_component, PlotComponent):
                ylims = self._set_y_limits(sub_component.x_values,
                                           sub_component.y_values,
                                           margin=margin)
                file_index = 0
                if color_by == "file":
                    file_index = comp_index*len(sub_component.labels)

                for index, (label, y_value, y_error)\
                        in enumerate(zip(sub_component.labels,
                                         sub_component.y_values,
                                         sub_component.y_errors)):
                    color = colors[file_index+index]
                    if self.readouts[0].band == "lband":
                        ax_left, ax_right = plot_broken_axis(
                            ax, sub_component.x_values,
                            y_value, y_error, self.lband_bounds,
                            self.mband_bounds, ax_left, ax_right,
                            color=color, ylims=ylims, label=label,
                            error=error, no_fill=no_fill,
                            err_percentile=offset*comp_index)

                        d = .015
                        kwargs_diagonal = dict(transform=ax_left.transAxes,
                                               color="k", clip_on=False)
                        ax_left.plot((1-d, 1+d), (-d, +d), **kwargs_diagonal)
                        ax_left.plot((1-d, 1+d), (1-d, 1+d), **kwargs_diagonal)

                        kwargs_diagonal.update(transform=ax_right.transAxes)
                        ax_right.plot((-d, +d), (1-d, 1+d), **kwargs_diagonal)
                        ax_right.plot((-d, +d), (-d, +d), **kwargs_diagonal)

                        ax.set_xlabel(xlabel, labelpad=20)
                        ax_left.set_ylabel(name)

                        handles.append(mlines.Line2D(
                            [], [], color=color, label=label))

                    else:
                        ax.plot(sub_component.x_values, y_value,
                                label=label, color=color)
                        if error:
                            if no_fill:
                                no_fill_index = int(
                                    np.ceil(sub_component.x_values.shape[-1]*comp_index*offset))
                                ax.errorbar(
                                    sub_component.x_values[no_fill_index],
                                    y_value[no_fill_index],
                                    yerr=np.mean(y_error),
                                    color=color, capsize=3)
                            else:
                                ax.fill_between(sub_component.x_values,
                                                y_value+y_error, y_value-y_error,
                                                color=color, alpha=0.2)

                        ax.set_ylim(*ylims)
                        ax.set_xlabel(xlabel)
                        ax.set_ylabel(name)

                if show_legend:
                    if handles:
                        if "left" in legend_location:
                            ax_legend = ax_left
                        else:
                            ax_legend = ax_right
                        ax_legend.legend(fontsize=legend_size, handles=handles,
                                         loc=legend_location, framealpha=0.5)
                    else:
                        ax.legend(fontsize=legend_size,
                                  loc=legend_location, framealpha=0.5)
            else:
                sub_component(ax, color_by=color_by, **kwargs)

        kwargs_layout = {"pad": 3.0, "h_pad": 2.0, "w_pad": 4.0}\
            if self.readouts[0].band == "lband" else {}
        plt.tight_layout(**kwargs_layout)
        return ax_left

    # TODO: Sharex, sharey and should be added
    def plot(self, save: Optional[bool] = False,
             subplots: Optional[bool] = False,
             figsize: Optional[List[int]] = None,
             dpi: Optional[int] = 100,
             sharex: Optional[bool] = False,
             rax: Optional[bool] = False,
             **kwargs) -> Optional[Tuple[Figure, Axes]]:
        """Combines the individual components into one plot.

        The size and dimension of the plot is automatically determined

        Parameters
        ----------
        save: bool, optional
            If toggled, saves the plot to the self.save_path file with the
            self.plot_name
        subplots : bool, optional
            If True, the components will be plotted in subplots.
        sharex : bool, optional
            If True, the x-axis will be shared.
        error : bool, optional
            If True, the error bars will be plotted.
        rax : bool, optional
            If True, then the figure and axes will be returned.
        kwargs: dict, optional

        Returns
        -------
        matplotlib.figure.Figure, optional
            The figure of the plot.
        matplotlib.axes.Axes, optional
            The axes of the plot.
        """
        if subplots:
            columns = self.num_components
            rows = len(list(self.components.values())[0]["values"])
        else:
            columns = 1 if self.num_components == 1 else\
                (3 if self.num_components >= 3 else 2)
            rows = np.ceil(self.num_components/columns).astype(int)\
                if self.num_components != 1 else 1

        to_px = 1/dpi
        size = figsize if figsize is not None else OPTIONS.plot.size
        fig, axarr = plt.subplots(rows, columns, tight_layout=True,
                                  figsize=(size*to_px*columns, size*to_px*rows))

        if self.num_components != 1:
            if subplots:
                for comp_index, (name, component) in enumerate(self.components.items()):
                    for index, sub_component in enumerate(component["values"]):
                        if isinstance(sub_component, PlotComponent):
                            sub_component = [sub_component]
                        else:
                            sub_component = [partial(sub_component, readouts=[self.readouts[index]])]
                        axarr[index, comp_index] = self.plot_components(
                            axarr[index, comp_index], name, sub_component,
                            **component["kwargs"], **kwargs)
            else:
                for index, (ax, (name, component)) in enumerate(zip(
                        axarr.flatten(), self.components.items())):
                    axarr[index] = self.plot_components(
                        ax, name, component["values"],
                        **component["kwargs"], **kwargs)
        else:
            name, component = map(
                lambda x: x[0], zip(*self.components.items()))
            axarr = self.plot_components(
                axarr, name, component["values"], **component["kwargs"], **kwargs)

        if rax:
            return fig, axarr

        if save:
            plt.savefig(self.save_dir / self.plot_name,
                        format=Path(self.plot_name).suffix[1:])
        else:
            plt.show()
        plt.close()
