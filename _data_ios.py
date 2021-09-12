import numpy as np
import astropy.io.ascii as asc

from utils import SRWUtils


class SRWSubIOS:
    """
    This subclass contains functions designed for reading datasets from files and writing them back.
    """

    def __init__(self):
        self.pars = None
        self.catalog = None
        self.raw_flux = None
        self.raw_magn = None
        self.raw_merr = None
        self.clr_magn = None
        self.clr_merr = None
        self.vs_mask = None

    def read_catalog_file(self):
        self.catalog = asc.read(self.pars["catalog_file"], format="commented_header", delimiter="\t",
                                fill_values=[(asc.masked, "nan")])

    def write_catalog_file(self):
        asc.write(self.catalog, self.pars["catalog_file"], overwrite=True, delimiter="\t",
                  format="commented_header", fill_values=[(asc.masked, "nan")])

    def read_raw_files(self):
        images_number = len(SRWUtils.get_image_list(self.pars["image_dir"], self.pars["image_filter"]))
        stars_number = len(self.catalog)
        with open(self.pars["raw_flux_file"], "r") as flux_file,\
                open(self.pars["raw_magn_file"], "r") as magn_file,\
                open(self.pars["raw_merr_file"], "r") as merr_file:
            self.raw_flux = np.zeros((images_number, stars_number))
            self.raw_magn = np.zeros((images_number, stars_number))
            self.raw_merr = np.zeros((images_number, stars_number))
            for image_index in range(images_number):
                self.raw_flux[image_index] = np.array(flux_file.readline().split(sep="\t"))
                self.raw_magn[image_index] = np.array(magn_file.readline().split(sep="\t"))
                self.raw_merr[image_index] = np.array(merr_file.readline().split(sep="\t"))

    def write_raw_files(self):
        with open(self.pars["raw_flux_file"], "w") as flux_file, \
                open(self.pars["raw_magn_file"], "w") as magn_file, \
                open(self.pars["raw_merr_file"], "w") as merr_file:
            for image_index in range(self.raw_flux.shape[0]):
                flux_file.write("\t".join(self.raw_flux[image_index].astype(str)) + "\n")
                magn_file.write("\t".join(self.raw_magn[image_index].astype(str)) + "\n")
                merr_file.write("\t".join(self.raw_merr[image_index].astype(str)) + "\n")

    def read_clear_files(self):
        images_number = len(SRWUtils.get_image_list(self.pars["image_dir"], self.pars["image_filter"]))
        stars_number = len(self.catalog)
        with open(self.pars["clr_magn_file"], "r") as magn_file,\
                open(self.pars["clr_merr_file"], "r") as merr_file:
            self.clr_magn = np.zeros((images_number, stars_number))
            self.clr_merr = np.zeros((images_number, stars_number))
            for image_index in range(images_number):
                self.clr_magn[image_index] = np.array(magn_file.readline().split(sep="\t"))
                self.clr_merr[image_index] = np.array(merr_file.readline().split(sep="\t"))

    def write_clear_files(self):
        with open(self.pars["clr_magn_file"], "w") as magn_file,\
                open(self.pars["clr_merr_file"], "w") as merr_file:
            for image_index in range(self.clr_magn.shape[0]):
                magn_file.write("\t".join(self.clr_magn[image_index].astype(str)) + "\n")
                merr_file.write("\t".join(self.clr_merr[image_index].astype(str)) + "\n")

    def read_vs_mask(self):
        with open(self.pars["vs_mask_file"], "r") as mask_file:
            self.vs_mask = np.array(mask_file.readline().split("\t"), dtype=bool)

    def write_vs_mask(self):
        with open(self.pars["vs_mask_file"], "w") as mask_file:
            mask_file.write("\t".join(self.vs_mask.astype(str)))
