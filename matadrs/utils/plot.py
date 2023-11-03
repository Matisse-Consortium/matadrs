from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Tuple, List, Union, Optional

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation
from matplotlib.axes import Axes

from .readout import ReadoutFits
from .tools import unwrap_phases, calculate_uv_points
from .options import OPTIONS


# TODO: Add proper docs
def make_uv_tracks(ax, uv_coord: np.ndarray[float],
                   baselines: List[np.ndarray],
                   sta_label: List[np.ndarray],
                   declination: float, flag: bool,
                   symbol: str, color: str, sel_wl: float,
                   airmass_lim: float, show_text: bool) -> None:
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
    u_coords, v_coords = uv_coord
    latitude_paranal = EarthLocation.of_site(
        "paranal").geodetic.lat.to(u.rad)
    hamax = np.arccos(abs((1./airmass_lim-np.sin(latitude_paranal)
                           * np.sin(declination))/(np.cos(latitude_paranal)
                                                   * np.cos(declination))))
    hour_angles = np.linspace(-hamax, hamax, 1000)
    u_coord_tracks, v_coord_tracks = calculate_uv_points(baselines,
                                                         hour_angles,
                                                         latitude_paranal,
                                                         declination)

    ax.plot(u_coord_tracks, v_coord_tracks, '-', color='grey', alpha=0.5)
    ax.plot(-u_coord_tracks, -v_coord_tracks, '-', color='grey', alpha=0.5)
    ax.plot([0.], [0.], '+k', markersize=5, markeredgewidth=2, alpha=0.5)

    ax.plot(u_coords, v_coords, symbol, color=color,
            markersize=10, markeredgewidth=3)
    ax.plot(-u_coords, -v_coords, symbol,
            color=color, markersize=10, markeredgewidth=3)
    if show_text:
        ax.text(-u_coords-3.5, -v_coords-1.5, sta_label,
                fontsize="small", color='0', alpha=0.8)


# TODO: Rewrite all of this to encompass errors and suplots as well
@dataclass
class PlotComponent:
    """Class containing the elements required for a plot
    from the """
    labels: List = None
    x_values: List = None
    y_values: List = None
    y_errors: List = None


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

    def __init__(self, fits_files: Union[Path, List[Path]],
                 flux_files: Optional[Union[Path, List[Path]]] = None,
                 plot_name: Optional[str] = None,
                 save_path: Optional[Path] = None) -> None:
        """The class's constructor"""
        self.fits_files = [fits_files]\
            if not isinstance(fits_files, List) else fits_files

        if flux_files is None:
            self.flux_files = [None]*len(self.fits_files)
        elif isinstance(flux_files, List)\
                and len(flux_files) != len(self.fits_files):
            raise IOError("Flux files must either be None or be "
                          "the same number the (.fits)-files provided")
        else:
            self.flux_files = flux_files

        if save_path is None:
            self.save_path = Path("").cwd()
        else:
            self.save_path = Path(save_path)

        if plot_name is not None:
            self.plot_name = plot_name
        elif isinstance(fits_files, List) and len(fits_files) > 1:
            self.plot_name = "combined_fits"
        else:
            self.plot_name = f"{Path(fits_files).stem}"

        self.readouts = [ReadoutFits(fits_file, flux_file)
                         for fits_file, flux_file
                         in zip(self.fits_files, self.flux_files)]
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
                      data: List[np.ndarray], margin: Optional[float] = 0.05) -> Tuple[int, int]:
        """Sets the y-limits from the data with some margin"""
        if np.min(wavelength) >= 6:
            indices = np.where((wavelength > 8.5) | (wavelength < 12.5))
        else:
            indices_low = np.where((wavelength <= 4.8) & (wavelength >= 4.5))
            indices_high = np.where((wavelength >= 3.) & (wavelength <= 3.8))
            indices = np.hstack((indices_low, indices_high))
        ymin, ymax = data[:, indices].min(), data[:, indices].max()
        spacing = np.linalg.norm(ymax-ymin)*margin
        return ymin-spacing, ymax+spacing

    def plot_uv(self, ax: Axes, symbol: Optional[str] = "x",
                sel_wl: Optional[float] = None,
                airmass_lim: Optional[float] = 2.,
                show_text: Optional[List] = False) -> None:
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
        uv_max = 0
        for index, readout in enumerate(self.readouts):
            uv_coords = readout.oi_vis["UVCOORD"]
            if uv_max < (tmp_uv_max := uv_coords.max()):
                uv_max = tmp_uv_max
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
                make_uv_tracks(ax, uv_coord, baselines[uv_index],
                               sta_labels[uv_index], readout.dec*np.pi/180,
                               flags[uv_index], symbol,
                               OPTIONS["plot.colors"][index],
                               sel_wl, airmass_lim, show_text)

            xlabel, ylabel = "$u$ (m)", "$v$ (m)"
            uv_extent = int(uv_max + uv_max*0.25)

            ax.set_xlim([uv_extent, -uv_extent])
            ax.set_ylim([-uv_extent, uv_extent])
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

    def make_component(self, data_name: str,
                       legend_format: Optional[str] = "long",
                       unwrap: Optional[bool] = False,
                       period: Optional[int] = 360) -> Union[Callable, PlotComponent]:
        """Generates a pandas DataFrame that has all the plots' information

        Parameters
        ----------
        data_name : str
            The name of the data to be plotted. Determines the legend- and plot
            labels.
        legend_format : str, optional
            Sets the format of the legend: For all information set "long" and
            for only station names set "short".
        unwrap : bool, optional
        period : int, optional

        Returns
        -------
        component : list of either Callable or PlotComponent
        """
        component = []
        for readout in self.readouts:
            sub_component = PlotComponent(x_values=readout.oi_wl["EFF_WAVE"].data.squeeze())
            if data_name == "flux":
                sub_component.y_values = readout.oi_flux["FLUXDATA"]
                sub_component.y_errors = readout.oi_flux["FLUXERR"]
                sub_component.labels = readout.oi_array["TEL_NAME"]\
                    if len(sub_component.y_values) > 1 else ["Averaged"]
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
                sub_component.labels = labels
                if data_name == "vis":
                    sub_component.y_values = readout.oi_vis["VISAMP"]
                    sub_component.y_errors = readout.oi_vis["VISAMPERR"]
                elif data_name == "diff":
                    diff_phases = readout.oi_vis["VISPHI"]
                    diff_phases_err = readout.oi_vis["VISPHIERR"]
                    if unwrap:
                        diff_phases, diff_phases_err = unwrap_phases(diff_phases,
                                                                     diff_phases_err, period)
                    sub_component.y_values = diff_phases
                    sub_component.y_errors = diff_phases_err
                elif data_name == "vis2":
                    sub_component.y_values = readout.oi_vis2["VIS2DATA"]
                    sub_component.y_errors = readout.oi_vis2["VIS2ERR"]
                else:
                    raise KeyError("No data-type of that data name exists!")
            elif data_name == "cphases":
                cphases, cphases_err = readout.oi_t3["T3PHI"], readout.oi_t3["T3PHIERR"]
                if unwrap:
                    cphases, cphases_err = unwrap_phases(cphases,
                                                         cphases_err, period)
                sub_component.labels = readout.oi_t3["TRIANGLE"]
                sub_component.y_values = cphases
                sub_component.y_errors = cphases_err
            elif data_name == "uv":
                sub_component = self.plot_uv
            else:
                raise KeyError("Input data name cannot be queried!")
            component.append(sub_component)
        return component

    def add_flux(self, **kwargs):
        """Adds the total flux(es) as a subplot

        See Also
        --------
        """
        self.components["Total Flux [Jy]"] =\
            self.make_component("flux", **kwargs)
        return self

    def add_vis(self, corr_flux: Optional[bool] = False, **kwargs):
        """Adds the visibilities/correlated fluxes as a subplot"""
        label = "Correlated Flux [Jy]" if corr_flux else "Visibility"
        self.components[label] =\
            self.make_component("vis", **kwargs)
        return self

    def add_vis2(self, corr_flux: Optional[bool] = False, **kwargs):
        """Adds the squared visibilities as a subplot"""
        self.components["Squared visibilities"] =\
            self.make_component("vis2", **kwargs)
        return self

    def add_cphases(self, **kwargs):
        """Adds the closure phases as a subplot"""
        self.components[r"Closure phases [$^{\circ}$]"] =\
            self.make_component("cphases", **kwargs)
        return self

    def add_diff_phases(self, **kwargs):
        """Adds the differential phases as a subplot"""
        self.components[r"Differential phases [$^{\circ}$]"] =\
            self.make_component("diff", **kwargs)
        return self

    def add_uv(self, **kwargs):
        """Adds the (u, v)-coordinates as a subplot"""
        self.components["$(u, v)$-coordinates"] =\
            self.make_component("uv", **kwargs)
        return self

    def add_mosaic(self, **kwargs):
        """Combines multiple subplots to produce a mosaic plot"""
        self.add_uv(**kwargs).add_vis(corr_flux=True, **kwargs).add_vis2(**kwargs)
        self.add_flux(**kwargs).add_cphases(**kwargs).add_diff_phases(**kwargs)
        return self

    def plot_component(self, ax, name: str,
                       component: Union[Callable, PlotComponent],
                       no_xlabel: Optional[bool] = False,
                       error: Optional[bool] = False,
                       margin: Optional[float] = 0.05,
                       **kwargs) -> None:
        """Plots all the data of a single component.

        Parameters
        ----------
        ax :
        name : str
        component : callable or plotcomponent
        no_xlabel : bool, optional
        error : bool, optional
        margin : bool, optional
        kwargs: dict
        """
        xlabel = r"$\lambda$ [$\mathrm{\mu}$m]" if not no_xlabel else ""
        for index, sub_component in enumerate(component):
            if isinstance(sub_component, PlotComponent):
                for label, y_value, y_error in zip(sub_component.labels,
                                                   sub_component.y_values,
                                                   sub_component.y_errors):
                    ax.plot(sub_component.x_values, y_value, label=label)
                    if error:
                        ax.fill_between(sub_component.x_values,
                                        y_value+y_error, y_value-y_error, alpha=0.2)
                    ax.legend(fontsize="xx-small", loc="upper right", framealpha=0.5)
                test = self._set_y_limits(sub_component.x_values,
                                          sub_component.y_values,
                                          margin=margin)
                ax.set_ylim(*test)
                ax.set_xlabel(xlabel)
                ax.set_ylabel(name)
            else:
                sub_component(ax, **kwargs)

    # TODO: Sharex, sharey and subplots should be added
    def plot(self, save: Optional[bool] = False,
             subplots: Optional[bool] = False,
             sharex: Optional[bool] = False,
             format: Optional[str] = "pdf", **kwargs):
        """Combines the individual components into one plot.

        The size and dimension of the plot is automatically determined

        Parameters
        ----------
        save: bool, optional
            If toggled, saves the plot to the self.save_path file with the
            self.plot_name
        subplots : bool, optional
        sharex : bool, optional
        kwargs: dict, optional
        """
        columns = 1 if self.num_components == 1 else\
            (3 if self.num_components >= 3 else 2)
        rows = np.ceil(self.num_components/columns).astype(int)\
            if self.num_components != 1 else 1
        to_px = 1/plt.rcParams["figure.dpi"]
        fig, axarr = plt.subplots(rows, columns,
                                  constrained_layout=True,
                                  figsize=(512*to_px*columns, 512*to_px*rows))

        if self.num_components != 1:
            for ax, (name, component) in zip(axarr.flatten(), self.components.items()):
                self.plot_component(ax, name, component, **kwargs)
        else:
            name, component = map(
                lambda x: x[0], zip(*self.components.items()))
            self.plot_component(axarr, name, component, **kwargs)
        fig.tight_layout()

        if save:
            plt.savefig(self.save_path / (self.plot_name + f".{format}"),
                        format=format)
        else:
            plt.show()
        plt.close()
