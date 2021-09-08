from datetime import datetime

import numpy as np
import photutils as pht
from astropy.wcs import WCS
from astropy.stats import SigmaClip

from utils import SRWUtils


class SRWSubAPP:
    """
    This subclass contains functions designed for conducting aperture photometry.
    """

    def __init__(self):
        self.pars = None
        self.catalog = None
        self.raw_flux = None
        self.raw_magn = None
        self.raw_merr = None

    def aperture_photometry(self):
        image_list = SRWUtils.get_image_list(self.pars["image_dir"], self.pars["image_filter"])
        images_number = len(image_list)
        stars_number = len(self.catalog)

        self.raw_flux = np.zeros((images_number, stars_number))
        self.raw_magn = np.zeros((images_number, stars_number))
        self.raw_merr = np.zeros((images_number, stars_number))

        for image_index in range(images_number):
            image_header, image_data = SRWUtils.read_fits_file(self.pars["image_dir"] + image_list[image_index])
            wcs_object = WCS(image_header)

            stars_xy_coords = wcs_object.wcs_world2pix(self.catalog["RAJ2000"], self.catalog["DEJ2000"], 0)
            bad_stars_x_mask = np.where((stars_xy_coords[0] < self.pars["edge"] / 10) |
                                        (stars_xy_coords[0] > (image_data.shape[1] - self.pars["edge"] / 10)))[0]
            bad_stars_y_mask = np.where((stars_xy_coords[1] < self.pars["edge"] / 10) |
                                        (stars_xy_coords[1] > (image_data.shape[0] - self.pars["edge"] / 10)))[0]
            bad_stars_mask = np.concatenate((bad_stars_x_mask, bad_stars_y_mask), axis=0)
            if len(bad_stars_mask) > 0:
                stars_xy_coords[0][bad_stars_mask] = 0
                stars_xy_coords[1][bad_stars_mask] = 0

            background_object = pht.Background2D(image_data, (100, 100),
                                                 filter_size=(10, 10),
                                                 sigma_clip=SigmaClip(sigma=3.),
                                                 bkg_estimator=pht.MedianBackground())
            image_data = image_data - background_object.background
            signal_sky = background_object.background_rms_median

            x_centroids_coord, y_centroids_coord = \
                pht.centroids.centroid_sources(image_data,
                                               xpos=stars_xy_coords[0],
                                               ypos=stars_xy_coords[1],
                                               box_size=self.pars["search_box"])
            apertures_object = pht.CircularAperture(np.vstack((x_centroids_coord,
                                                               y_centroids_coord)).T,
                                                    r=self.pars["aperture"])

            flux_image_data = np.array(pht.aperture_photometry(image_data,
                                                               apertures_object,
                                                               method="exact")["aperture_sum"])
            magn_image_data = image_header["ZEROPOI"] - 2.5 * np.log10(flux_image_data) + 2.5 * np.log10(
                image_header["EXPTIME"])
            merr_image_data = 1.0857 * np.sqrt(flux_image_data * self.pars["gain"] + apertures_object.area * (
                    signal_sky * self.pars["gain"] + self.pars["ron"] ** 2)) / (flux_image_data * self.pars["gain"])

            flux_image_data[bad_stars_mask] = np.nan
            magn_image_data[bad_stars_mask] = np.nan
            merr_image_data[bad_stars_mask] = np.nan

            self.raw_flux[image_index] = flux_image_data
            self.raw_magn[image_index] = magn_image_data
            self.raw_merr[image_index] = merr_image_data

    @staticmethod
    def _logs_apphot(image_index, image_header):
        with open("LOGS.txt", "a") as logs_file:
            logs_file.write(f"{datetime.now()}\tAPPHOT: Processed image #{image_index + 1}\t"
                            f"DATE-OBS {image_header['DATE-OBS']}")
            print(f"{datetime.now()}\tAPPHOT: Processed image # {image_index + 1}\t"
                  f"DATE-OBS {image_header['DATE-OBS']}")
