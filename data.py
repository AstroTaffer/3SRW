from datetime import datetime


class SRWData:
    """
    This class is the one intended for usage.
    It inherits all methods from internal subclasses.
    New datasets are created as the corresponding functions are called.
    Datasets include:
        - Set of processing parameters
        - Internal catalog
        - Raw flux, magnitude and magnitude error
        - Clear flux, magnitude and magnitude error
        - [mask]

    Class dictionary pars_names contains full parameters names and can be used as a handbook.
    """

    pars_names = {"image_dir": "Images directory",
                  "raw_flux_file": "Raw flux file",
                  "raw_magn_file": "Raw magnitude file",
                  "raw_merr_file": "Raw magnitude error file",
                  "clr_magn_file": "Clear magnitude file",
                  "clr_merr_file": "Clear magnitude error file",
                  "catalog_file": "Internal catalog file",
                  "vss_mask_file": "Variable stars search mask file",
                  "image_filter": "Images filter",
                  # Only files with image_filter in their names are processed

                  "ext_catalog": "External catalog",
                  "mag_lim": "Maximum magnitude limit",
                  # Only stars brighter than mag_lim are included
                  "image_edge": "Images edge",

                  "aperture": "Aperture radius",
                  "search_box": "Centroid search box",
                  "gain": "CCD gain",
                  "ron": "CCD readout noise",

                  "mmd": "Ensemble stars maximal magnitude difference",
                  "isr": "Ensemble stars initial search radius in arcminutes",
                  "msr": "Ensemble stars maximum search radius in arcminutes",

                  "std_lim": "Ensemble stars standart deviation limit",

                  "vss_method": "Variable stars search method"}

    def __init__(self, image_dir, **kwargs):
        # noinspection PyDictCreation
        self.pars = {}
        self.catalog = None
        self.raw_flux = None
        self.raw_magn = None
        self.raw_merr = None
        self.vs_mask = None

        # These parameters are used in output files names definition
        # And therefore must be specified out of turn
        self.pars["image_filter"] = kwargs.get("image_filter", "V")
        self.pars["aperture"] = kwargs.get("aperture", 4)
        # To me: set aperture to 4, 6 or 8
        self.pars["vss_method"] = kwargs.get("vss_method", "ALL")
        # Set to "ALL", "SIGMA", "ERRORS" or "ROMS2"

        self.pars["image_dir"] = image_dir
        self.pars["catalog_file"] = kwargs.get("catalog_file", "CATALOG.txt")
        self.pars["raw_flux_file"] = kwargs.get("raw_flux_file", f"RAW_FLUX_\
{self.pars['image_filter']}{self.pars['aperture']}.txt")
        self.pars["raw_magn_file"] = kwargs.get("raw_magn_file", f"RAW_MAGN_\
{self.pars['image_filter']}{self.pars['aperture']}.txt")
        self.pars["raw_merr_file"] = kwargs.get("raw_merr_file", f"RAW_MERR_\
{self.pars['image_filter']}{self.pars['aperture']}.txt")
        self.pars["clr_magn_file"] = kwargs.get("clr_magn_file", f"CLR_MAGN_\
{self.pars['image_filter']}{self.pars['aperture']}.txt")
        self.pars["clr_merr_file"] = kwargs.get("clr_merr_file", f"CLR_MERR_\
{self.pars['image_filter']}{self.pars['aperture']}.txt")
        self.pars["vss_mask_file"] = kwargs.get("vss_mask_file", f"VSS_MASK_\
{self.pars['image_filter']}{self.pars['aperture']}_{self.pars['vss_method']}.txt")

        self.pars["ext_catalog"] = kwargs.get("ext_catalog", "II/336")
        self.pars["mag_lim"] = kwargs.get("mag_lim", 14)
        self.pars["image_edge"] = kwargs.get("image_edge", 100)

        self.pars["search_box"] = kwargs.get("search_box", 15)
        self.pars["gain"] = kwargs.get("gain", 1.3)
        self.pars["ron"] = kwargs.get("ron", 10)

        self.pars["mmd"] = kwargs.get("mmd", 2.0)
        self.pars["isr"] = kwargs.get("isr", 5)
        self.pars["msr"] = kwargs.get("msr", 30)

        self.pars["std_lim"] = kwargs.get("std_lim", 3)

        self._logs_pars()

    def _logs_pars(self):
        with open("LOGS.txt", "a") as logs_file:
            logs_file.write(2 * "\n" + 20 * "-" + "\n")
            logs_file.write(f"Current date and time\t{datetime.now()}\n")
            logs_file.write("New set of parameters is loaded:\n")
            print(2 * "\n" + 20 * "-" + "\n")
            print(f"Current date and time\t{datetime.now()}\n")
            print("New set of parameters is loaded:\n")

            for key in self.pars.keys():
                logs_file.write(f"{key}\t{SRWData.pars_names[key]} {self.pars[key]}\n")
                print(f"{key}\t{SRWData.pars_names[key]} {self.pars[key]}\n")
