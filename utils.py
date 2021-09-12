from os import listdir

import astropy.io.fits as fits


class SRWUtils:
    """
    This class contains multi purpose functions used in various places of this module.
    """

    @classmethod
    def get_image_list(cls, image_dir, image_filter):
        dir_content = listdir(image_dir)
        image_list = []

        for file in dir_content:
            if file.count(".fits") and file.count(image_filter):
                image_list.append(file)

        cls._image_list_logs(image_list, dir_content)

        return image_list

    @staticmethod
    def _image_list_logs(image_list, dir_content):
        with open("LOGS.txt", "a") as logs_file:
            logs_file.write(f"Image directory scanned\nSelected {len(image_list)} images out of {len(dir_content)}\n")
            print(f"Image directory scanned\nSelected {len(image_list)} images out of {len(dir_content)}")

    @staticmethod
    def read_fits_file(file_path):
        with fits.open(file_path) as hdu_list:
            hdu_list.verify("fix")
            image_header = hdu_list[0].header
            image_data = hdu_list[0].data

        return image_header, image_data
