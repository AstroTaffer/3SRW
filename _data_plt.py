from datetime import datetime

import numpy as np
import astropy.coordinates as coord
import astropy.units as units
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta

from utils import SRWUtils


class SRWSubMPL:
    """
    This subclass contains functions designed for plotting various graphs.
    """

    def __init__(self):
        self.pars = None
        self.catalog = None
        self.clr_magn = None
        self.clr_merr = None
        self.vs_mask = None

    def plot_vs_light_curves(self):
        image_list = SRWUtils.get_image_list(self.pars["image_dir"], self.pars["image_filter"])
        images_number = len(image_list)
        stars_number = len(self.catalog)
        bjd_list = np.zeros(images_number)

        for _ in range(images_number):
            header = SRWUtils.read_fits_file(self.pars["image_dir"] + image_list[_])
            obs_time = Time(header["DATE-OBS"], format="fits") + TimeDelta(header["EXPTIME"], format="sec") / 2
            # noinspection PyUnresolvedReferences
            sc_obj = coord.SkyCoord(header["CRVAL1"],
                                    header["CRVAL2"],
                                    unit=(units.deg, units.deg),
                                    frame="icrs")
            # noinspection PyUnresolvedReferences
            el_object = coord.EarthLocation.from_geodetic(lon=header["LONGITUD"] * units.deg,
                                                          lat=header["LATITUDE"] * units.deg,
                                                          height=header["ALTITUDE"] * units.m)
            bjd_list[_] = (obs_time + obs_time.light_travel_time(sc_obj, location=el_object)).jd

        vs_counter = 0
        plt_fig, plt_ax = plt.subplots(dpi=150)
        xlabel_bjd = (bjd_list[0] // 1)
        bjd_list %= 1

        for _ in range(stars_number):
            if self.vs_mask[_]:
                plt_ax.set_xlabel(f"BJD - {xlabel_bjd}")
                plt_ax.set_ylabel("m")
                plt_ax.grid()
                plt_ax.plot(bjd_list, self.clr_magn[:, _], "k.", markersize=2)
                plt.savefig(f"LC_{self.pars['image_filter']}{self.pars['aperture']}_"
                            f"{self.pars['vss_method']}_{_ + 1}.png")
                plt_ax.cla()
                vs_counter += 1

            if (_ + 1) % 100 == 0:
                self._logs_light_curves(vs_counter)

        self._logs_light_curves(vs_counter)

    def plot_merr_graph(self):
        mean_star_magn = np.nanmean(self.clr_magn, axis=0)
        mean_star_merr = np.nanmean(self.clr_merr, axis=0)

        plt_fig, plt_ax = plt.subplots(dpi=150)
        plt_ax.set_xlabel(r"$\langle m \rangle$")
        plt_ax.set_ylabel(r"$\langle m_{err} \rangle$")
        plt_ax.grid()
        plt_ax.plot(np.nanmean(self.clr_magn, axis=0), np.nanmean(self.clr_merr, axis=0), "k.", markersize=2)
        plt.savefig(f"MERR(MAGN)_{self.pars['image_filter']}{self.pars['aperture']}.png")

        for _ in range(len(self.catalog)):
            if self.vs_mask[_]:
                plt_ax.plot(mean_star_magn[_], mean_star_merr[_], "bo", markersize=4)

        if self.pars["vss_method"] == "ERRORS":
            bad_stars_indexes = np.where(np.isnan(mean_star_magn) | np.isnan(mean_star_merr))
            mean_star_magn = np.delete(mean_star_magn, bad_stars_indexes)
            mean_star_merr = np.delete(mean_star_merr, bad_stars_indexes)

            apr_merr_func = np.polyfit(mean_star_magn, np.log(mean_star_merr), deg=1)
            # Cut w=np.sqrt(merr_mean) - might be harmful
            apr_xpoints = np.linspace(np.min(mean_star_magn), np.max(mean_star_magn), num=100)

            plt_ax.plot(apr_xpoints, np.exp(apr_merr_func[0] * apr_xpoints + apr_merr_func[1]),
                        "r-", linewidth=1)
            plt_ax.plot(apr_xpoints, np.exp(apr_merr_func[0] * apr_xpoints + apr_merr_func[1]) + 0.0005,
                        "g-", linewidth=1)
            # 0.0005 is the width of the curve. It was determined empirically.
            # That's not good - need more accurate and general method

            plt.savefig(f"MERR(MAGN)_{self.pars['image_filter']}{self.pars['aperture']}_ERRORS.png")

        elif self.pars["vss_method"] == "ROMS2":
            pass

    @staticmethod
    def _logs_light_curves(vs_counter):
        with open("LOGS.txt", "a") as logs_file:
            logs_file.write(f"{datetime.now()}\tVSPLT: Plotted {vs_counter} light curves\n")
            print(f"{datetime.now()}\tVSPLT: Plotted {vs_counter} light curves\n")
