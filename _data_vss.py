from datetime import datetime

import numpy as np

from utils import SRWUtils


class SRWSubVSS:
    """
    This subclass contains functions designed for variable stars search.
    """

    def __init__(self):
        self.pars = None
        self.catalog = None
        self.vs_mask = None
        self.clr_magn = None
        self.clr_merr = None

    def vs_search(self):
        if self.pars["vss_method"] == "ALL":
            self._vs_all()
        elif self.pars["vss_method"] == "SIGMA":
            self._vs_sigma()
        elif self.pars["vss_method"] == "ERRORS":
            self._vs_merr()
        elif self.pars["vss_method"] == "ROMS2":
            self._vs_roms2()

    def _vs_all(self):
        stars_number = len(self.catalog)
        self.vs_mask = np.zeros(stars_number, dtype=bool)
        vs_counter = 0

        for _ in range(stars_number):
            if np.isfinite(self.clr_magn[:, _]).all():
                self.vs_mask[_] = True
                vs_counter += 1

            if (_ + 1) % 100 == 0:
                self._logs_vss(_, vs_counter)

        self._logs_vss(stars_number, vs_counter)

    def _vs_sigma(self):
        stars_number = len(self.catalog)
        self.vs_mask = np.zeros(stars_number, dtype=bool)
        vs_counter = 0

        mean_star_magn = np.nanmean(self.clr_magn, axis=0)
        std_star_magn = np.nanstd(self.clr_magn, axis=0)

        finite_id = np.isfinite(mean_star_magn)
        apr_std_func = np.polyfit(mean_star_magn[finite_id], std_star_magn[finite_id], deg=1)

        for _ in range(stars_number):
            # Use width
            if std_star_magn[_] - (apr_std_func[0] * mean_star_magn[_] + apr_std_func[1]) > 0.20:
                self.vs_mask[_] = True
                vs_counter += 1
                # 0.20 is the width of the curve. It was determined empirically.
                # That's not good - need more accurate and general method

            if (_ + 1) % 100 == 0:
                self._logs_vss(_, vs_counter)

        self._logs_vss(stars_number, vs_counter)

    def _vs_merr(self):
        stars_number = len(self.catalog)
        self.vs_mask = np.zeros(stars_number, dtype=bool)
        vs_counter = 0

        mean_star_magn = np.nanmean(self.clr_magn, axis=0)
        mean_star_merr = np.nanmean(self.clr_merr, axis=0)

        finite_id = np.isfinite(mean_star_magn)
        apr_merr_func = np.polyfit(mean_star_magn[finite_id], np.log(mean_star_merr[finite_id]),
                                   deg=1, w=np.sqrt(mean_star_merr))

        for _ in range(stars_number):
            if mean_star_merr[_] - np.exp(apr_merr_func[0] * mean_star_magn[_] + apr_merr_func[1]) > 0.0005:
                self.vs_mask[_] = True
                vs_counter += 1
                # 0.0005 is the width of the curve. It was determined empirically.
                # That's not good - need more accurate and general method

            if (_ + 1) % 100 == 0:
                self._logs_vss(_, vs_counter)

        self._logs_vss(stars_number, vs_counter)

    def _vs_roms2(self):
        stars_number = len(self.catalog)
        images_number = len(SRWUtils.get_image_list(self.pars["image_dir"], self.pars["image_filter"]))
        self.vs_mask = np.zeros(stars_number, dtype=bool)
        vs_counter = 0

        mean_star_magn = np.nanmean(self.clr_magn, axis=0)
        std_star_magn = np.nanstd(self.clr_magn, axis=0)

        good_id = np.where(mean_star_magn > 10)[0]
        apr_std_func = np.polyfit(mean_star_magn[good_id], np.log(std_star_magn[good_id]),
                                  deg=1)
        # Cut w=np.sqrt(std_star_magn[good_id]) - harmful

        for _ in range(stars_number):
            med_star_magn = np.nanmedian(self.clr_magn[:, _])
            etha = 0.0
            for __ in range(images_number):
                etha += abs((self.clr_magn[__, _] - med_star_magn) /
                            (apr_std_func[0] * self.clr_magn[__, _] + apr_std_func[1]))
            etha /= (images_number - 1)

            if etha > 1:
                self.vs_mask[_] = True
                vs_counter += 1

            if (_ + 1) % 100 == 0:
                self._logs_vss(_, vs_counter)

        self._logs_vss(stars_number, vs_counter)

    @staticmethod
    def _logs_vss(all_stars, var_stars):
        with open("LOGS.txt", "a") as logs_file:
            logs_file.write(f"{datetime.now()}\tVSS: Processed {all_stars} stars, "
                            f"{var_stars} recognised as variable\n")
            print(f"{datetime.now()}\tVSS: Processed {all_stars} stars, {var_stars} recognised as variable")
