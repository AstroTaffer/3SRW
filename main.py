from data import SRWData


image_path = "H:/CODING_LAIR/IMAGES/EAST_V/"
daob = SRWData(image_path, aperture=6, image_filter="V")

# daob.create_catalog()
# daob.write_catalog_file()
daob.read_catalog_file()

daob.aperture_photometry()
daob.write_raw_files()
# daob.read_raw_files()

daob.differential_photometry()
daob.write_clear_files()
# daob.read_clear_files()


# ALL
daob.pars["vss_method"] = "ALL"
daob.vs_search()
daob.write_vs_mask()
daob.plot_vs_light_curves()


"""
# ERRORS
daob.pars["vss_method"] = "ERRORS"
daob.vs_search()
daob.write_vs_mask()
daob.plot_merr_graph()
daob.plot_vs_light_curves()
"""

"""
# SIGMA
daob.pars["vss_method"] = "SIGMA"
daob.vs_search()
daob.write_vs_mask()
daob.plot_sigma_graph()
daob.plot_vs_light_curves()
"""

"""
# ROMS2
daob.pars["vss_method"] = "ROMS2"
daob.vs_search()
daob.write_vs_mask()
daob.plot_ln_sigma_graph()
daob.plot_sigma_graph()
daob.plot_vs_light_curves()
"""
