"""$Id: $

This file is part of the Matisse pipeline GUI series
Copyright (C) 2017- Observatoire de la Côte d'Azur

Created on Sat Mar 17 06:39:49 2018
@author: fmillour
fmillour@oca.eu

This software is a computer program whose purpose is to show oifits
files from the MATISSE instrument.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

Changelog:
2018-03-23: New functions: oi_data_select_frame, filter_oi_list, open_oi_dir, show_vis2_tf2_vs_time, show_oi_vs_time (jvarga)
2018-03-26: New GUI interface ready: oi_data_select_frame (jvarga)
2018-04-04: Updated GUI and extended functionality: input file/folder textbox, filter for target name,
            more bands (JHK) available (for e.g. AMBER data), plot with or without errorbars, plot V or V2 (jvarga)
2023-12-06: Updated saving of plots (automated folder creation) and documentation of functions (mScheuck)
"""
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
from matplotlib import pyplot as plt

from ..utils import robust


def show_seeing(
        list_of_dicts: List[Dict],
        saveplots: Optional[bool] = False,
        output_path: Optional[Path] = None,
        show: Optional[bool] = True, **kwargs) -> None:
    """Plots the seeing for the oifits-files.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    saveplots : bool, optional
        If True, the plots are saved in the output_path. The default is False.
    output_path : pathlib.Path, optional
        Path to save the plots. The default is None.
    show : bool, optional
        If True, the plots are shown. The default is True.
    """
    if list_of_dicts:
        seeing_arr = []
        tau0_arr = []
        mjd_arr = []

        for dic in list_of_dicts:
            try:
                mjd_arr.append(np.array(dic["VIS2"]["TIME"])[0])
            except KeyError:
                continue
            seeing_arr.append(dic["SEEING"])
            tau0_arr.append(dic["TAU0"] * 1.0e3)
        fig = plt.figure(figsize=(7, 7))
        axs2 = fig.add_subplot(2, 1, 2)
        axs1 = fig.add_subplot(2, 1, 1, sharex=axs2)
        plt.setp(axs1.get_xaxis().get_offset_text(), visible=False)
        axs1.tick_params(axis="x", which="both", bottom=False,
                         top=False, labelbottom=False)

        axs2.plot(mjd_arr, tau0_arr, "s")
        print(f"Tau0_MEAN = {np.mean(tau0_arr)} +- {np.std(tau0_arr)}")
        axs2.set_ylabel("Tau0 (ms)")
        axs2.set_xlabel(r"$\mathrm{MJD}$")

        axs1.plot(mjd_arr, seeing_arr, "d")
        print(f"Seeing_MEAN = {np.mean(seeing_arr)} +- {np.std(seeing_arr)}")
        axs1.set_ylabel("Seeing")

        axs1.set_ylim(0.0, np.max(seeing_arr) + 0.5)
        axs2.set_ylim(0.0, np.max(tau0_arr) + 1.0)

        for dic in list_of_dicts:
            target_names = dic["TARGET"]
            try:
                axs1.text(
                    np.array(dic["VIS2"]["TIME"])[0],
                    np.max(seeing_arr) + 0.5 + 0.1,
                    target_names.replace("_", " "),
                    rotation=90, va="bottom", fontsize=8)
            except KeyError:
                continue

        plt.subplots_adjust(hspace=0)

        if saveplots:
            label = "_SEEING"
            fig.savefig(output_path / f"{label}.png", dpi=150)
            fig.savefig(output_path / f"{label}.eps", format="eps", dpi=300)
            plt.close(fig)
        if show:
            plt.show()


def show_vis_vs_seeing(
        list_of_dicts: List[Dict], wlenRange: List[float],
        numPeak: Optional[int] = 1,
        numPeak2: Optional[int] = 6,
        saveplots: Optional[bool] = False,
        vis: Optional[bool] = True,
        output_dir: Optional[Path] = None,
        show: Optional[bool] = True) -> None:
    """Plots the visibility vs. seeing for the oifits-files.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    wlenRange : list
        Range of wavelength to plot.
    numPeak : int, optional
        Number of peaks to plot. The default is 1.
    numPeak2 : int, optional
        Number of peaks to plot. The default is 6.
    saveplots : bool, optional
        If True, the plots are saved in the output_path. The default is False.
    output_path : pathlib.Path, optional
        Path to save the plots. The default is None.
    show : bool, optional
        If True, the plots are shown. The default is True.
    """
    if list_of_dicts:
        target_names_TF2 = []
        TF2_BCD_arr = []
        TF2_MJD_arr = []
        TF2_arr = []
        TF2err_arr = []
        TF_arr = []
        TFerr_arr = []
        TF2_arr_2nd = []
        TF2err_arr_2nd = []
        TF_arr_2nd = []
        TFerr_arr_2nd = []
        TF2_sta_index = []
        TF2_seeing = []
        TF2_tau0 = []
        TF2_pwv = []

        for dic in list_of_dicts:
            wl = np.array(dic["WLEN"])
            wlenRange_idx = np.logical_and(
                wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6
            )
            if sum(wlenRange_idx) > 0:
                category = dic["CATEGORY"].lower()
                if "cal" in category:
                    try:
                        datay_vis = np.array(dic["TF2"]["TF"])
                        datayerr_vis = np.array(dic["TF2"]["TFERR"])
                        datay = np.array(dic["TF2"]["TF2"])
                        datayerr = np.array(dic["TF2"]["TF2ERR"])
                        datax = np.array(dic["TF2"]["TIME"])
                        n_rows = datay.shape[0]

                        for i in range(n_rows // 6):
                            if dic["BCD1NAME"] == "IN":
                                BCD1 = 1
                            elif dic["BCD1NAME"] == "OUT":
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic["BCD2NAME"] == "IN":
                                BCD2 = 1
                            elif dic["BCD2NAME"] == "OUT":
                                BCD2 = 0
                            else:
                                BCD2 = 1
                            TF2_arr.append(
                                robust.mean(datay[i * 6 + numPeak - 1, wlenRange_idx])
                            )
                            TF2err_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            TF_arr.append(
                                np.sqrt(3.4)
                                * robust.mean(
                                    datay_vis[i * 6 + numPeak - 1, wlenRange_idx]
                                )
                            )
                            TFerr_arr.append(
                                np.sqrt(3.4)
                                * robust.mean(datayerr_vis[i, wlenRange_idx])
                            )
                            TF2_arr_2nd.append(
                                robust.mean(datay[i * 6 + numPeak2 - 1, wlenRange_idx])
                            )
                            TF2err_arr_2nd.append(
                                robust.mean(datayerr[i, wlenRange_idx])
                            )
                            TF_arr_2nd.append(
                                np.sqrt(3.4)
                                * robust.mean(
                                    datay_vis[i * 6 + numPeak2 - 1, wlenRange_idx]
                                )
                            )
                            TFerr_arr_2nd.append(
                                np.sqrt(3.4)
                                * robust.mean(datayerr_vis[i, wlenRange_idx])
                            )
                            TF2_BCD_arr.append([BCD1, BCD2])
                            TF2_MJD_arr.append(datax[i])
                            TF2_seeing.append(dic["SEEING"])
                            TF2_tau0.append(dic["TAU0"])
                            TF2_pwv.append(dic["PWV"])
                            target_names_TF2.append(dic["TARGET"])
                            sta_index = np.sort(dic["TF2"]["STA_INDEX"][i])
                            TF2_sta_index.append(sta_index)
                    except:
                        print(dic["TARGET"], dic["DATEOBS"], "No CAL TF2 data found.")

        target_names_TF2 = np.array(target_names_TF2)
        TF2_BCD_arr = np.array(TF2_BCD_arr)
        print(TF2_BCD_arr)
        ind_BCD = np.where(TF2_BCD_arr[:, 0] == 0)
        print("ind_BCD = {0}".format(ind_BCD))
        TF2_MJD_arr = np.array(TF2_MJD_arr)
        TF2_arr = np.array(TF2_arr)
        TF2err_arr = np.array(TF2err_arr)
        TF_arr = np.array(TF_arr)
        TFerr_arr = np.array(TFerr_arr)
        TF2_arr_2nd = np.array(TF2_arr_2nd)
        TF2err_arr_2nd = np.array(TF2err_arr_2nd)
        TF_arr_2nd = np.array(TF_arr_2nd)
        TFerr_arr_2nd = np.array(TFerr_arr_2nd)
        TF2_sta_index = np.array(TF2_sta_index)
        TF2_seeing = np.array(TF2_seeing)
        TF2_tau0 = np.array(TF2_tau0)
        TF2_pwv = np.array(TF2_pwv)

        ind_tau04 = np.where(TF2_tau0 * 1e3 < 4)
        ind_tau08 = np.where((TF2_tau0 * 1e3 > 4) & (TF2_tau0 * 1e3 < 8))
        ind_tau0sup = np.where(TF2_tau0 * 1e3 > 8)

        ind_seeing06 = np.where(TF2_seeing < 0.6)
        ind_seeing1 = np.where((TF2_seeing > 0.6) & (TF2_seeing < 1.0))
        ind_seeingsup = np.where(TF2_seeing > 1.0)

        ind_tau0 = np.argsort(TF2_tau0)
        TF2_tau0_sorted = TF2_tau0[ind_tau0]
        TF2_arr_sorted = TF2_arr[ind_tau0]
        TF_arr_sorted = TF_arr[ind_tau0]

        fig = plt.figure(figsize=(20, 20))
        axs2 = fig.add_subplot(2, 1, 1)
        axs1 = fig.add_subplot(2, 1, 2)

        if vis:
            axs2.plot(
                TF2_tau0[ind_seeingsup] * 1e3,
                TF_arr[ind_seeingsup],
                "o",
                label='seeing > 1"',
            )
            axs2.plot(
                TF2_tau0[ind_seeing1] * 1e3,
                TF_arr[ind_seeing1],
                "o",
                label='0.6" < seeing < 1"',
            )
            axs2.plot(
                TF2_tau0[ind_seeing06] * 1e3,
                TF_arr[ind_seeing06],
                "o",
                label='seeing < 0.6"',
            )
        else:
            axs2.plot(
                TF2_tau0[ind_seeingsup] * 1e3,
                np.sqrt(TF2_arr[ind_seeingsup]),
                "o",
                markersize=12,
                label='seeing > 1"',
            )
            axs2.plot(
                TF2_tau0[ind_seeing1] * 1e3,
                np.sqrt(TF2_arr[ind_seeing1]),
                "o",
                markersize=12,
                label='0.6" < seeing < 1"',
            )
            axs2.plot(
                TF2_tau0[ind_seeing06] * 1e3,
                np.sqrt(TF2_arr[ind_seeing06]),
                "o",
                markersize=12,
                label='seeing < 0.6"',
            )

        axs2.set_xlabel("$\\tau_0$ [ms]", fontsize=36)
        axs2.set_ylabel("V", fontsize=36)

        if vis:
            axs1.plot(
                TF2_seeing[ind_tau04],
                TF_arr[ind_tau04],
                "o",
                label=r"2ms < $\tau_0$ < 4ms",
            )
            axs1.plot(
                TF2_seeing[ind_tau08],
                TF_arr[ind_tau08],
                "o",
                label=r"4ms < $\tau_0$ < 8ms",
            )
            axs1.plot(
                TF2_seeing[ind_tau0sup],
                TF_arr[ind_tau0sup],
                "o",
                label=r"$\tau_0$ > 8ms",
            )
        else:
            axs1.plot(
                TF2_seeing[ind_tau04],
                np.sqrt(TF2_arr[ind_tau04]),
                "o",
                markersize=12,
                label=r"$\tau_0$ < 4ms",
            )
            axs1.plot(
                TF2_seeing[ind_tau08],
                np.sqrt(TF2_arr[ind_tau08]),
                "o",
                markersize=12,
                label=r"4ms < $\tau_0$ < 8ms",
            )
            axs1.plot(
                TF2_seeing[ind_tau0sup],
                np.sqrt(TF2_arr[ind_tau0sup]),
                "o",
                markersize=12,
                label=r"$\tau_0$ > 8ms",
            )
        axs1.set_ylabel("V", fontsize=36)
        axs1.set_xlabel('Seeing ["]', fontsize=36)

        axs1.set_ylim(0, 1)
        axs2.set_xlim(1, 14)
        axs2.set_ylim(0, 1)
        axs1.legend(loc=0, fontsize=30)
        axs2.legend(loc=4, fontsize=30)
        axs1.tick_params(axis="both", labelsize="36")
        axs2.tick_params(axis="both", labelsize="36")
        plt.tight_layout()
        fig.subplots_adjust(hspace=0.2)

        if saveplots:
            label = "_V_WRT_SEEING"
            fig.savefig(output_dir / f"{label}.png", dpi=300)
            fig.savefig(output_dir / f"{label}.eps", format="eps", dpi=300)
            plt.close(fig)
        if show:
            plt.show()


def show_vis_vs_seeing_GRA4MAT(
        list_of_dicts: List[Dict],
        wlenRange: List[float],
        numPeak: int = 1,
        saveplots: Optional[bool] = False,
        vis: Optional[bool] = True,
        output_path: Optional[Path] = None,
        show: Optional[bool] = True) -> None:
    """Plots the vis for the oifits-files vs seeing.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    wlenRange : list of float
        The wavelength range to plot.
    numPeak : int, optional
        The number of peaks to plot. The default is 1.
    saveplots : bool, optional
        If True, the plots are saved in the output_path. The default is False.
    vis : bool, optional
        If True, the plots are saved in the output_path. The default is True.
    output_path : pathlib.Path, optional
        Path to save the plots. The default is None.
    show : bool, optional
        If True, the plots are shown. The default is True.
    """
    if list_of_dicts:
        target_names_TF2 = []
        TF2_BCD_arr = []
        TF2_MJD_arr = []

        TF2_arr = []
        TF2err_arr = []
        TF_arr = []
        TFerr_arr = []
        TF2_arr_gra = []
        TF2err_arr_gra = []
        TF_arr_gra = []
        TFerr_arr_gra = []
        TF2_sta_index = []
        TF2_seeing = []
        TF2_tau0 = []
        TF2_sta_index_gra = []
        TF2_seeing_gra = []
        TF2_tau0_gra = []

        for dic in list_of_dicts:
            wl = np.array(dic["WLEN"])
            wlenRange_idx = np.logical_and(
                wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6
            )
            if sum(wlenRange_idx) > 0:
                category = dic["CATEGORY"].lower()
                if "cal" in category:
                    try:
                        datay_vis = np.array(dic["TF2"]["TF"])
                        datayerr_vis = np.array(dic["TF2"]["TFERR"])
                        datay = np.array(dic["TF2"]["TF2"])
                        datayerr = np.array(dic["TF2"]["TF2ERR"])
                        datax = np.array(dic["TF2"]["TIME"])
                        n_rows = datay.shape[0]

                        if dic["FT"] == "MATISSE":
                            for i in range(n_rows // 6):
                                if dic["BCD1NAME"] == "IN":
                                    BCD1 = 1
                                elif dic["BCD1NAME"] == "OUT":
                                    BCD1 = 0
                                else:
                                    BCD1 = 0
                                if dic["BCD2NAME"] == "IN":
                                    BCD2 = 1
                                elif dic["BCD2NAME"] == "OUT":
                                    BCD2 = 0
                                else:
                                    BCD2 = 1
                                TF2_arr.append(
                                    robust.mean(
                                        datay[i * 6 + numPeak - 1, wlenRange_idx]
                                    )
                                )
                                TF2err_arr.append(
                                    robust.mean(datayerr[i, wlenRange_idx])
                                )
                                TF_arr.append(
                                    np.sqrt(3.4)
                                    * robust.mean(
                                        datay_vis[i * 6 + numPeak - 1, wlenRange_idx]
                                    )
                                )
                                TFerr_arr.append(
                                    np.sqrt(3.4)
                                    * robust.mean(datayerr_vis[i, wlenRange_idx])
                                )
                                TF2_seeing.append(dic["SEEING"])
                                TF2_tau0.append(dic["TAU0"])
                                target_names_TF2.append(dic["TARGET"])
                                sta_index = np.sort(dic["TF2"]["STA_INDEX"][i])
                                TF2_sta_index.append(sta_index)

                        elif dic["FT"] == "GRAVITY":
                            for i in range(n_rows // 6):
                                if dic["BCD1NAME"] == "IN":
                                    BCD1 = 1
                                elif dic["BCD1NAME"] == "OUT":
                                    BCD1 = 0
                                else:
                                    BCD1 = 0
                                if dic["BCD2NAME"] == "IN":
                                    BCD2 = 1
                                elif dic["BCD2NAME"] == "OUT":
                                    BCD2 = 0
                                else:
                                    BCD2 = 1
                                TF2_arr_gra.append(
                                    robust.mean(
                                        datay[i * 6 + numPeak - 1, wlenRange_idx]
                                    )
                                )
                                TF2err_arr_gra.append(
                                    robust.mean(datayerr[i, wlenRange_idx])
                                )
                                TF_arr_gra.append(
                                    np.sqrt(3.4)
                                    * robust.mean(
                                        datay_vis[i * 6 + numPeak - 1, wlenRange_idx]
                                    )
                                )
                                TFerr_arr_gra.append(
                                    np.sqrt(3.4)
                                    * robust.mean(datayerr_vis[i, wlenRange_idx])
                                )
                                TF2_seeing_gra.append(dic["SEEING"])
                                TF2_tau0_gra.append(dic["TAU0"])
                                target_names_TF2_gra.append(dic["TARGET"])
                                sta_index = np.sort(dic["TF2"]["STA_INDEX"][i])
                                TF2_sta_index_gra.append(sta_index)
                    except Exception:
                        print(dic["TARGET"], dic["DATEOBS"], "No CAL TF2 data found.")

        TF2_arr = np.array(TF2_arr)
        TF2err_arr = np.array(TF2err_arr)
        TF_arr = np.array(TF_arr)
        TFerr_arr = np.array(TFerr_arr)
        TF2_arr_gra = np.array(TF2_arr_gra)
        TF2err_arr_gra = np.array(TF2err_arr_gra)
        TF_arr_gra = np.array(TF_arr_gra)
        TFerr_arr_gra = np.array(TFerr_arr_gra)

        TF2_sta_index = np.array(TF2_sta_index)
        TF2_seeing = np.array(TF2_seeing)
        TF2_tau0 = np.array(TF2_tau0)
        TF2_sta_index_gra = np.array(TF2_sta_index_gra)
        TF2_seeing_gra = np.array(TF2_seeing_gra)
        TF2_tau0_gra = np.array(TF2_tau0_gra)

        fig = plt.figure(figsize=(20, 20))
        axs2 = fig.add_subplot(2, 1, 1)
        axs1 = fig.add_subplot(2, 1, 2)

        if vis:
            axs1.plot(
                TF2_seeing, TF_arr, "o",
                markersize=12, label="MATISSE standalone"
            )
            axs1.plot(TF2_seeing_gra, TF_arr_gra,
                      "o", markersize=12, label="GRA4MAT")
        else:
            axs1.plot(
                TF2_seeing,
                np.sqrt(TF2_arr),
                "o",
                markersize=12,
                label="MATISSE standalone",
            )
            axs1.plot(
                TF2_seeing_gra,
                np.sqrt(TF2_arr_gra),
                "o",
                markersize=12,
                label="GRA4MAT",
            )
        axs1.set_xlabel('Seeing ["]', fontsize=36)
        axs1.set_ylabel("V", fontsize=36)

        if vis:
            axs2.plot(
                TF2_tau0 * 1e3, TF_arr, "o", markersize=12, label="MATISSE standalone"
            )
            axs2.plot(
                TF2_tau0_gra * 1e3, TF_arr_gra, "o", markersize=12, label="GRA4MAT"
            )
        else:
            axs2.plot(
                TF2_tau0 * 1e3,
                np.sqrt(TF2_arr),
                "o",
                markersize=12,
                label="MATISSE standalone",
            )
            axs2.plot(
                TF2_tau0_gra * 1e3,
                np.sqrt(TF2_arr_gra),
                "o",
                markersize=12,
                label="GRA4MAT",
            )
        axs2.set_xlabel("$\\tau_0$ [ms]", fontsize=36)
        axs2.set_ylabel("V", fontsize=36)

        axs1.set_ylim(0, 1)
        axs2.set_xlim(1, 10)
        axs2.set_ylim(0, 1)

        axs1.legend(loc=0, fontsize=30)
        axs2.legend(loc=4, fontsize=30)

        axs1.tick_params(axis="both", labelsize="36")
        axs2.tick_params(axis="both", labelsize="36")
        plt.tight_layout()
        fig.subplots_adjust(hspace=0.2)

        if saveplots:
            label = "_V_WRT_SEEING"
            fig.savefig(output_path / f"{label}.png", dpi=300)
            fig.savefig(output_path / f"{label}.eps", format="eps", dpi=300)
            plt.close(fig)
        if show:
            plt.show()


###############################################################################
def show_clo_vs_seeing(
        list_of_dicts: List[Dict],
        wlenRange: List[float],
        numPeak: Optional[int] = 1,
        saveplots: Optional[bool] = False,
        output_path: Optional[Path] = None,
        show: Optional[bool] = True) -> None:
    """Plots the closure vs. seeing for the oifits-files.

    Parameters
    ----------
    list_of_dicts : list of dict
        List of oifits-files.
    wlenRange : list
        Range of wavelength to plot.
    numPeak : int, optional
        Number of peaks to plot. The default is 1.
    saveplots : bool, optional
        If True, the plots are saved in the output_path. The default is False.
    output_path : pathlib.Path, optional
        Path to save the plots. The default is None.
    show : bool, optional
        If True, the plots are shown. The default is True.
    """
    if list_of_dicts:
        target_names_PHI3 = []
        PHI3_BCD_arr = []
        PHI3_MJD_arr = []
        PHI3_arr = []
        PHI3err_arr = []
        PHI3_sta_index = []
        PHI3_seeing = []
        PHI3_tau0 = []

        for dic in list_of_dicts:
            wl = np.array(dic["WLEN"])
            wlenRange_idx = np.logical_and(
                wl > wlenRange[0] / 1.0e6, wl < wlenRange[1] / 1.0e6
            )
            if sum(wlenRange_idx) > 0:
                category = dic["CATEGORY"].lower()
                if "cal" in category:
                    try:
                        datay = np.array(dic["T3"]["CLOS"])
                        datayerr = np.array(dic["T3"]["CLOSERR"])
                        datax = np.array(dic["T3"]["TIME"])
                        n_rows = datay.shape[0]

                        for i in range(n_rows // 4):
                            if dic["BCD1NAME"] == "IN":
                                BCD1 = 1
                            elif dic["BCD1NAME"] == "OUT":
                                BCD1 = 0
                            else:
                                BCD1 = 0
                            if dic["BCD2NAME"] == "IN":
                                BCD2 = 1
                            elif dic["BCD2NAME"] == "OUT":
                                BCD2 = 0
                            else:
                                BCD2 = 0
                            PHI3_arr.append(
                                robust.mean(datay[i * 4 + numPeak - 1, wlenRange_idx])
                            )
                            PHI3err_arr.append(robust.mean(datayerr[i, wlenRange_idx]))
                            PHI3_BCD_arr.append([BCD1, BCD2])
                            PHI3_MJD_arr.append(datax[i])
                            PHI3_seeing.append(dic["SEEING"])
                            PHI3_tau0.append(dic["TAU0"])
                            target_names_PHI3.append(dic["TARGET"])
                            sta_index = np.sort(dic["T3"]["STA_INDEX"][i])
                            PHI3_sta_index.append(sta_index)
                    except Exception:
                        print(dic["TARGET"], dic["DATEOBS"], "No CAL TF2 data found.")

        target_names_PHI3 = np.array(target_names_PHI3)
        PHI3_BCD_arr = np.array(PHI3_BCD_arr)
        ind_BCD = np.where(PHI3_BCD_arr[:, 0] == 0)
        print(f"ind_BCD = {ind_BCD}")

        PHI3_MJD_arr = np.array(PHI3_MJD_arr)
        PHI3_arr = np.array(PHI3_arr)
        PHI3err_arr = np.array(PHI3err_arr)
        PHI3_sta_index = np.array(PHI3_sta_index)
        PHI3_seeing = np.array(PHI3_seeing)
        PHI3_tau0 = np.array(PHI3_tau0)

        PHI3_arr = PHI3_arr[ind_BCD]
        PHI3_tau0 = PHI3_tau0[ind_BCD]
        PHI3_seeing = PHI3_seeing[ind_BCD]

        ind_tau02 = np.where(PHI3_tau0 * 1e3 < 2)
        ind_tau04 = np.where((PHI3_tau0 * 1e3 > 2) & (PHI3_tau0 * 1e3 < 4))
        ind_tau06 = np.where((PHI3_tau0 * 1e3 > 4) & (PHI3_tau0 * 1e3 < 6))
        ind_tau08 = np.where((PHI3_tau0 * 1e3 > 6) & (PHI3_tau0 * 1e3 < 8))
        ind_tau0sup = np.where(PHI3_tau0 * 1e3 > 8)

        ind_seeing04 = np.where(PHI3_seeing < 0.4)
        ind_seeing06 = np.where((PHI3_seeing > 0.4) & (PHI3_seeing < 0.6))
        ind_seeing08 = np.where((PHI3_seeing > 0.6) & (PHI3_seeing < 0.8))
        ind_seeing1 = np.where((PHI3_seeing > 0.8) & (PHI3_seeing < 1.0))
        ind_seeingsup = np.where(PHI3_seeing > 1.0)

        print(f"median Seeing = {np.median(PHI3_seeing)}")
        print(f"median Tau0 = {np.median(PHI3_tau0)}")

        fig = plt.figure(figsize=(9, 9))
        axs3 = fig.add_subplot(3, 1, 3)
        axs2 = fig.add_subplot(3, 1, 2)
        axs1 = fig.add_subplot(3, 1, 1)

        axs3.plot(PHI3_seeing, PHI3_tau0 * 1e3, "d")
        axs3.set_xlabel("Seeing (as)")
        axs3.set_ylabel("Tau0 (ms)")

        # axs2.plot(TF2_tau0[ind_tau0]*1e+3,np.sqrt(TF2_arr[ind_tau0]),'s')
        # axs2.plot(TF2_tau0[ind_tau0]*1e+3,p(TF2_tau0[ind_tau0]*1e+3))
        axs2.plot(
            PHI3_tau0[ind_seeing04] * 1e3,
            PHI3_arr[ind_seeing04],
            "d",
            label="seeing < 0.4",
        )
        axs2.plot(
            PHI3_tau0[ind_seeing06] * 1e3,
            PHI3_arr[ind_seeing06],
            "d",
            label="0.4 < seeing < 0.6",
        )
        axs2.plot(
            PHI3_tau0[ind_seeing08] * 1e3,
            PHI3_arr[ind_seeing08],
            "d",
            label="0.6 < seeing < 0.8",
        )
        axs2.plot(
            PHI3_tau0[ind_seeing1] * 1e3,
            PHI3_arr[ind_seeing1],
            "d",
            label="0.8 < seeing < 1.0",
        )
        axs2.plot(
            PHI3_tau0[ind_seeingsup] * 1e3,
            PHI3_arr[ind_seeingsup],
            "d",
            label="seeing > 1.0",
        )

        axs2.set_xlabel("Tau0 (ms)")
        axs2.set_ylabel("Closure phase [deg]")

        axs1.plot(PHI3_seeing[ind_tau02], PHI3_arr[ind_tau02],
                  "d", label="tau0 < 2ms")
        axs1.plot(
            PHI3_seeing[ind_tau04], PHI3_arr[ind_tau04],
            "d", label="2ms < tau0 < 4ms"
        )
        axs1.plot(
            PHI3_seeing[ind_tau06], PHI3_arr[ind_tau06],
            "d", label="4ms < tau0 < 6ms"
        )
        axs1.plot(
            PHI3_seeing[ind_tau08], PHI3_arr[ind_tau08],
            "d", label="6ms < tau0 < 8ms"
        )
        axs1.plot(
            PHI3_seeing[ind_tau0sup], PHI3_arr[ind_tau0sup],
            "d", label="tau0 > 8ms"
        )
        axs1.set_ylabel("Closure phase [deg]")
        axs1.set_xlabel("Seeing (as)")
        axs1.set_title(
            f"BCD-calibrated Closure phase (bispectrum peak {numPeak})")

        axs1.set_ylim(-3, 3)
        axs2.set_ylim(-3, 3)
        axs1.legend(loc=0)
        axs2.legend(loc=4)
        plt.tick_params(axis="both", labelsize="12")
        plt.tight_layout()
        if saveplots:
            label = "_PHI3_WRT_SEEING"
            fig.savefig(output_path / f"{label}.png", dpi=150)
            fig.savefig(output_path / f"{label}.eps", format="eps", dpi=300)
            plt.close(fig)
        if show:
            plt.show()
