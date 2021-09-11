from datetime import datetime

import numpy as np


class SRWSubMPL:
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
        self.vs_mask = [False in range(stars_number)]
        vs_counter = 0

        for _ in range(stars_number):
            if not np.isnan(self.clr_magn[:, _]).all():
                self.vs_mask[_] = True
                vs_counter += 1

            if (_ + 1) % 100 == 0:
                self._logs_vss(_, vs_counter)

        self._logs_vss(stars_number, vs_counter)

    def _vs_sigma(self):
        pass

    def _vs_merr(self):
        stars_number = len(self.catalog)
        self.vs_mask = [False in range(stars_number)]
        vs_counter = 0

        mean_star_magn = np.nanmean(self.clr_magn, axis=0)
        mean_star_merr = np.nanmean(self.clr_merr, axis=0)
        bad_stars_indexes = np.where(np.isnan(mean_star_magn) | np.isnan(mean_star_merr))
        apr_mean_star_magn = np.delete(mean_star_magn, bad_stars_indexes)
        apr_mean_star_merr = np.delete(mean_star_merr, bad_stars_indexes)

        apr_merr_func = np.polyfit(apr_mean_star_magn, np.log(apr_mean_star_merr), deg=1)
        # Cut w=np.sqrt(merr_mean) - might be harmful

        for _ in range(stars_number):
            if mean_star_merr[_] - np.exp(apr_merr_func[0] * mean_star_magn[_] + apr_merr_func[1]) > 0.0005:
                self.vs_mask[_] = True
                vs_counter += 1

            if (_ + 1) % 100 == 0:
                self._logs_vss(_, vs_counter)

        self._logs_vss(stars_number, vs_counter)

    def _vs_roms2(self):
        pass

    @staticmethod
    def _logs_vss(all_stars, var_stars):
        with open("LOGS.txt", "a") as logs_file:
            logs_file.write(f"{datetime.now()}\tVSS: Processed {all_stars} stars, "
                            f"{var_stars} recognised as variable\n")
            print(f"{datetime.now()}\tDIFPHOT: Processed {all_stars} stars, {var_stars} recognised as variable\n")
