from datetime import datetime

import numpy as np
import astropy.coordinates as coord
import astropy.units as units

from utils import SRWUtils


class SRWSubDFP:
    """
    This subclass contains functions designed for conducting differential (ensemble) photometry.
    """

    def __init__(self):
        self.pars = None
        self.catalog = None
        self.raw_flux = None
        self.raw_magn = None
        self.raw_merr = None
        self.clr_magn = None
        self.clr_merr = None

    def differential_photometry(self):
        image_list = SRWUtils.get_image_list(self.pars["image_dir"], self.pars["image_filter"])
        images_number = len(image_list)
        stars_number = len(self.catalog)
        counter_success = 0

        self.clr_magn = np.zeros((images_number, stars_number))
        self.clr_merr = np.zeros((images_number, stars_number))

        # FIXME - if catalog is read from file then it needs to be multiplied by units.deg
        # FIXME - Otherwise, if calculated - does not
        # catalog_stars = coord.SkyCoord(ra=self.catalog["RAJ2000"], dec=self.catalog["DEJ2000"], frame="icrs")
        # noinspection PyUnresolvedReferences
        catalog_stars = coord.SkyCoord(ra=self.catalog["RAJ2000"] * units.deg,
                                       dec=self.catalog["DEJ2000"] * units.deg,
                                       frame="icrs")

        for target_star_index in range(stars_number):
            if np.isfinite(self.raw_magn[:, target_star_index]).any():

                ens_stars_indexes = []
                stars_goodness = [True for _ in range(stars_number)]
                search_radius = self.pars["isr"]
                # noinspection PyUnresolvedReferences
                target_star = coord.SkyCoord(ra=self.catalog[target_star_index]["RAJ2000"] * units.deg,
                                             dec=self.catalog[target_star_index]["DEJ2000"] * units.deg,
                                             frame="icrs")

                while len(ens_stars_indexes) < 10 and search_radius <= self.pars["msr"]:
                    # noinspection PyUnresolvedReferences
                    ens_stars_indexes = np.where(target_star.separation(catalog_stars) <
                                                 (search_radius * units.arcmin))[0]

                    check_index = 0
                    while check_index < len(ens_stars_indexes):
                        # if check_star is target_star or check_star is "all-NaN" or "not good"
                        # or has too different magnitude:
                        # delete check_star from ensemble
                        if ens_stars_indexes[check_index] == target_star_index or \
                                np.isnan(self.raw_magn[:, check_index]).all() or \
                                not stars_goodness[ens_stars_indexes[check_index]] or \
                                np.absolute(np.nanmean(self.raw_magn[:, target_star_index]) -
                                            np.nanmean(self.raw_magn[:, ens_stars_indexes[check_index]]) >
                                            self.pars["msr"]):
                            ens_stars_indexes = np.delete(ens_stars_indexes, check_index)
                        else:
                            check_index += 1

                    if len(ens_stars_indexes) < 10:
                        search_radius += 1
                        stars_goodness = [True for _ in range(stars_number)]
                        continue

                    ens_stars_magn = np.zeros((images_number, len(ens_stars_indexes)))
                    ens_stars_merr = np.zeros((images_number, len(ens_stars_indexes)))
                    for _ in range(len(ens_stars_indexes)):
                        ens_stars_magn[:, _] = self.raw_magn[:, ens_stars_indexes[_]]
                        ens_stars_merr[:, _] = self.raw_merr[:, ens_stars_indexes[_]]

                    eca_image_weight = np.ones(len(ens_stars_indexes)) / np.square(np.nanmean(ens_stars_merr,
                                                                                              axis=0))
                    eca_mean_image_magn = np.nansum(ens_stars_magn * eca_image_weight,
                                                    axis=1) / np.nansum(eca_image_weight)

                    eca_correction = eca_mean_image_magn - np.nanmean(eca_mean_image_magn)
                    ens_stars_magn -= eca_correction.reshape((-1, 1))

                    eca_std = np.nanstd(ens_stars_magn, axis=0)
                    if np.max(eca_std) > self.pars["std_lim"]:
                        stars_goodness[ens_stars_indexes[np.argmax(eca_std)]] = False
                        ens_stars_indexes = []
                        continue
                    else:
                        self.clr_magn[:, target_star_index] = self.raw_magn[:, target_star_index] - eca_correction
                        self.clr_merr[:, target_star_index] = np.sqrt(
                            (1 / np.nansum(1 / np.square(ens_stars_merr), axis=1))
                            + np.square(self.raw_merr[:, target_star_index]))
                        counter_success += 1

                if len(ens_stars_indexes) < 10:
                    self.clr_magn[:, target_star_index] = np.nan
                    self.clr_merr[:, target_star_index] = np.nan
            else:
                # It would be better to just outright delete 'all-NaN' stars so I don't have to worry about them later
                # But I don't have neither time nor will to worry about it right now

                self.clr_magn[:, target_star_index] = np.nan
                self.clr_merr[:, target_star_index] = np.nan
                continue

            if (target_star_index + 1) % 100 == 0:
                self._logs_diffphot(target_star_index + 1, counter_success)

        self._logs_diffphot(stars_number, counter_success)

    @staticmethod
    def _logs_diffphot(all_stars, success_stars):
    # FIXME Logs are somehow called in wrong way
        with open("LOGS.txt", "a") as logs_file:
            logs_file.write(f"{datetime.now()}\tDIFPHOT: Processed {all_stars} stars, "
                            f"{success_stars} successfully\n")
            print(f"{datetime.now()}\tDIFPHOT: Processed {all_stars} stars, {success_stars} successfully\n")
