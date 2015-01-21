# coding: utf-8

""" Script to run the filtering GUI and determine the optimal bandpass filter
    for the SkyMapper star. """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Third-party.
import numpy as np
from astropy.io import fits

# Module specific
import gui


# Load the SkyMapper star
filename = "/Users/arc/research/sm0313/uves/data/SM0313_RED_SCI_POINT_BLUE.fits"
image = fits.open(filename)
header, flux = image[0].header, image[0].data

dispersion = header["CRVAL1"] + np.arange(header["NAXIS1"]) * header["CDELT1"]

# Shouldn't have to do this for the SkyMapper data, but let's be sure:
if "LTV1" in header:
    dispersion -= header["LTV1"] * header["CDELT1"]
if header.get("CTYPE1", None) == "AWAV-LOG":
    dispersion = np.exp(dispersion)

# Run the interactive bandpass filter
bandpass_filter = gui.InteractiveFilter(dispersion, flux)
