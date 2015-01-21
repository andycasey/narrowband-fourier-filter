# coding: utf-8

""" Interactive GUI to visualise the effects of a narrow-band Fourier filter on
    astronomical spectra. """

from __future__ import division, print_function

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

__all__ = ["InteractiveFilter"]

# Standard library
import logging
import os

# Third-party
import numpy as np
import traitsui.api as ui
import traits.api as traits
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from wx import CallAfter

# Module-specific
from gui_core import (Figure, MaxNLocator, MPLFigureEditor)

def InteractiveFilter(dispersion, flux, *args, **kwargs):
    """ Return an instance of the interactive filter. """

    gui = _InteractiveFilter(dispersion, flux, *args, **kwargs)
    gui.configure_traits()
    return gui


class _InteractiveFilter(traits.HasTraits):

    # Smoothing.
    smoothing_kernel = traits.Float(0)
    apply_smoothing_button = traits.Button("Apply Smoothing")

    # Filtering.
    filtered = traits.Bool(False)
    bandwidth = traits.Float(100)
    apply_filter_button = traits.Button("Apply Filter")

    # Plotting.
    plot_spectra = traits.Instance(Figure)

    # Export.
    save_spectrum = traits.Button("Save")

    # Three figures, and a bandwidth slider at the bottom.
    view = ui.View(
        ui.VGroup(
            ui.Group(ui.Item("plot_spectra", editor=MPLFigureEditor(),
                show_label=False)),
            ui.HGroup(
                ui.Item("smoothing_kernel", label="Smoothing kernel (A)",
                    width=50, padding=5),
                ui.Item("apply_smoothing_button", show_label=False, padding=5,
                    enabled_when="smoothing_kernel > 0"),
                ui.spring,
                ui.Item("bandwidth", label="Bandwidth (A)", width=50, padding=5),
                ui.Item("apply_filter_button", show_label=False, padding=5),
                ui.spring,
                ui.Item("save_spectrum", show_label=False, padding=5,
                    enabled_when="filtered")
            )
        ), 
        width     = 1280,
        height    = 700,
        title     = "Fourier Filtering",
        resizable = True)

    
    def __init__(self, dispersion, flux, *args, **kwargs):

        dispersion, flux = map(np.array, (dispersion, flux))
        assert dispersion.size == flux.size

        self._dispersion = dispersion
        self._original_flux = flux
        self._smoothed_original_flux = flux.copy()
        self._filtered_flux = flux.copy()
        

    def _plot_spectra_default(self):
        """ Initialises the spectrum display. """

        figure = Figure()
        figure.subplots_adjust(hspace=0, wspace=0)

        # Plot the original spectra.
        self._axes = [figure.add_subplot(3, 1, 1)]

        # Plot the filtered spectra.
        self._axes.append(figure.add_subplot(3, 1, 2, sharex=self._axes[0],
            sharey=self._axes[0]))

        # Plot the filtered, normalised spectra.
        self._axes.append(figure.add_subplot(3, 1, 3, sharex=self._axes[0]))

        # Set y-axis major locations.
        self._axes[0].yaxis.set_major_locator(MaxNLocator(4))
        self._axes[1].yaxis.set_major_locator(MaxNLocator(4))
        self._axes[2].yaxis.set_major_locator(MaxNLocator(4))

        # Set axes labels.
        self._axes[0].set_ylabel("Original Flux")
        self._axes[1].set_ylabel("Filtered Flux")
        self._axes[2].set_ylabel("Normalised Flux")
        self._axes[2].set_xlabel("Wavelength")

        # Initialise the plot lines in each axes.
        self._axes[0].plot(self._dispersion, self._original_flux, c="k")
        # This will be to overplot the filtered spectrum.
        self._axes[0].plot([], [], c="r", lw=2, zorder=100)

        # Update the x-limits in the second and third plot. The first and second
        # plot have shared axes, so this will actually update them all.
        self._axes[0].set_xlim(self._dispersion[0], self._dispersion[-1])

        # Update the y-limits in the first and second axes. The third axis will
        # be normalised flux, so it will have different y-limits.
        self._axes[0].set_ylim(self._clip_ylimits(self._original_flux))
        
        # Here we do blank plots so that the object is accessible later on, but
        # the filtered and normalised spectra is actually drawn by other methods
        # that rely on the filter bandwidth.
        self._axes[1].plot([], [], c="r", lw=2)
        self._axes[2].plot([], [], c="k")
        self._axes[2].axhline(1.0, linestyle=":", c="#666666", zorder=-1)

        # Set the background.
        figure.patch.set_facecolor("w")
        figure.tight_layout()
        figure.subplots_adjust(hspace=0, wspace=0)

        return figure


    def _apply_smoothing_button_fired(self, _):
        """ Apply some level of smoothing to the data. """

        print("Smoothing the data..")

        # Smoothing kernel is in Angstroms. Need it in pixels.
        # px = Angstroms * (pixels/Angstrom)
        pixel_size = np.diff(self._dispersion).mean() # Angstroms/pixel
        kernel_size = self.smoothing_kernel / pixel_size
        self._smoothed_original_flux = gaussian_filter(self._original_flux,
            kernel_size)

        # Redraw the original and normalised spectra because they are smoothed.
        self._draw_original_spectrum()
        self._draw_normalised_spectrum()
        CallAfter(self.plot_spectra.canvas.draw)


    def _apply_filter_button_fired(self, _):
        """ Apply the new filter to the data. """

        print("Applying filter..")

        self.filtered = True
        self._apply_bandpass_filter()
        self._draw_filtered_spectrum()
        self._draw_normalised_spectrum()
        CallAfter(self.plot_spectra.canvas.draw)


    def _apply_bandpass_filter(self):
        """ Apply the bandpass filter to the original data. """

        print("Appling Fourier filter of bandwidth {0:.1f} Angstoms to the data"
            .format(self.bandwidth))

        size = self._original_flux.size
        _ = np.empty(3 * size)
        _[:size] = self._original_flux[::-1].copy()
        _[size:2*size] = self._original_flux.copy()
        _[2*size:] = self._original_flux[::-1].copy()

        # Bandwidth is in Angstroms, but we need it in pixels.
        pixel_size = np.diff(self._dispersion).mean() # Angstroms/pixel
        pixel_bandwidth = self.bandwidth / pixel_size

        fft = np.fft.rfft(_)
        x = np.arange(fft.size, dtype=float)
        model = np.fft.irfft(fft * x / (pixel_bandwidth + x))[size:2*size]
        self._filtered_flux = gaussian_filter(self._original_flux - model,
            size**0.5)


    def _draw_filtered_spectrum(self):
        """ Draw the filtered spectrum. """

        # Any smoothing?
        # [TODO]
        # Draw it on the top axes.
        self._axes[0].lines[1].set_data(self._dispersion, self._filtered_flux)

        # And draw it on it's own in the middle axis.
        self._axes[1].lines[0].set_data(self._dispersion, self._filtered_flux)


    def _draw_original_spectrum(self):
        """ Draw the original spectrum. """

        self._axes[0].lines[0].set_data(self._dispersion,
            self._smoothed_original_flux)


    def _draw_normalised_spectrum(self):
        """ Draw the normalised spectrum. """

        normalised_flux = self._smoothed_original_flux/self._filtered_flux
        self._axes[2].lines[0].set_data(self._dispersion, normalised_flux)

        # Set y-limits
        self._axes[2].set_ylim(0, 1.2)


    def _clip_ylimits(self, data):
        """ Some rule for default y-limits. """

        ideal = [np.nanmin(data), np.nanmedian(data) + 3 * np.nanstd(data)]
        return np.clip(ideal, 0, np.inf)


    def _save_spectrum_fired(self, _):
        """ Save the normalised spectrum to disk. """

        normalised_flux = self._original_flux/self._filtered_flux
        with open("normalised-spectrum.pickle", "w") as fp:
            pickle.dump((self._dispersion, normalised_flux), fp, -1)

        print("Saved unsmoothed, normalised data to normalised-spectrum.pickle")




if __name__ == "__main__":

    disp = np.arange(1000)
    flux = np.arange(1000) * 2

    # Load the SkyMapper star

    gui = InteractiveFilter(disp, flux)
    gui.configure_traits()


