import numpy as np
import astropy.coordinates as coord
import astropy.units as units
from astropy.wcs import WCS
from astroquery.vizier import Vizier

from utils import SRWUtils


class SRWSubCCT:
    """
    This subclass contains function designed for creating internal catalog.
    """

    def __init__(self):
        self.pars = None
        self.catalog = None

    def create_catalog(self):
        image = SRWUtils.get_image_list(self.pars["image_dir"], self.pars["image_filter"])[0]
        header = SRWUtils.read_fits_file(self.pars["image_dir"] + image)[0]

        wcs_object = WCS(header)
        image_center_crd = wcs_object.wcs_pix2world(header["NAXIS1"] / 2, header["NAXIS2"] / 2, 0)
        # Here only CD1_1 value was used because this simple method provides tolerable accuracy
        image_radius = abs(header["CD1_1"]) * ((header["NAXIS1"] ** 2 + header["NAXIS2"] ** 2) ** 0.5) / 2

        # noinspection PyUnresolvedReferences
        skycoord_object = coord.SkyCoord(ra=image_center_crd[0] * units.deg,
                                         dec=image_center_crd[1] * units.deg, frame="icrs")
        # noinspection PyUnresolvedReferences
        search_angle = coord.Angle(image_radius * units.deg)

        vizier_object = Vizier(columns=["RAJ2000", "DEJ2000", "Vmag"],
                               column_filters={"Vmag": "<" + self.pars["mag_lim"]})
        vizier_object.ROW_LIMIT = -1
        vizier_result = vizier_object.query_region(skycoord_object, radius=search_angle, catalog="II/336")
        self.catalog = vizier_result[0]

        stars_xy_coords = wcs_object.wcs_world2pix(self.catalog["RAJ2000"], self.catalog["DEJ2000"], 0)
        bad_stars_x_mask = np.where((stars_xy_coords[0] < self.pars["image_edge"]) |
                                    (stars_xy_coords[0] > (header["NAXIS1"] - self.pars["image_edge"])))[0]
        bad_stars_y_mask = np.where((stars_xy_coords[1] < self.pars["image_edge"]) |
                                    (stars_xy_coords[1] > (header["NAXIS2"] - self.pars["image_edge"])))[0]
        bad_stars_mask = np.concatenate((bad_stars_x_mask, bad_stars_y_mask), axis=0)
        if len(bad_stars_mask) > 0:
            self.catalog.remove_rows(bad_stars_mask)

        self.catalog.add_column(np.arange(1, len(self.catalog) + 1), index=0, name="ID")

        self._logs_catalog(image_center_crd, image_radius)

    def _logs_catalog(self, icc, ir):
        with open("LOGS.txt", "a") as logs_file:
            logs_file.write(self.catalog.info())
            logs_file.write(f"Image center: RA = {icc[0]}\tDEC = {icc[1]}\n")
            logs_file.write(f"Image radius: R = {ir}\n")
            print(self.catalog.info())
            print(f"Image center: RA = {icc[0]}\tDEC = {icc[1]}")
            print(f"Image radius: R = {ir}")
