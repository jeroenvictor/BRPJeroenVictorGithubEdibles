import numpy as np
from astropy.io import fits
import astropy.constants as cst
import matplotlib.pyplot as plt
from edibles.edibles_settings import datadir


class EdiblesSpectrum:
    # This object will contain a spectrum from EDIBLES,
    # and a set of methods to operate on the data.

    def loadSpectrum(self):
        # Assume the file is a DR3 product here.
        hdu = fits.open(self.filename)
        self.header = hdu[0].header
        self.target = self.header["OBJECT"]
        self.date = self.header["DATE-OBS"]
        self.flux = hdu[0].data
        self.flux_units = "arbitrary"
        crval1 = self.header["CRVAL1"]
        cdelt1 = self.header["CDELT1"]
        nwave = len(self.flux)
        grid = np.arange(0, nwave, 1)
        self.wave = (grid) * cdelt1 + crval1
        self.wave_units = "AA"
        self.reference_frame = "geocentric"
        self.v_bary = self.header["HIERARCH ESO QC VRAD BARYCOR"]
        self.bary_wave = self.wave + (self.v_bary / cst.c.to("km/s").value) * self.wave

    def __init__(self, filename):
        """
        Filename is relative to the DR3 directory
        """
        self.filename = datadir + filename
        self.loadSpectrum()

    def getSpectrum(self, xmin=None, xmax=None, bary=False):

        if (xmin is not None) and (xmax is not None):
            assert xmin < xmax, "xmin must be less than xmax"
            idx_tell = (self.wave > xmin) * (self.wave < xmax)
            if bary is True:
                idx_bary = (self.bary_wave > xmin) * (self.bary_wave < xmax)
                return self.bary_wave[np.where(idx_bary)], self.flux[np.where(idx_bary)]
            return self.wave[np.where(idx_tell)], self.flux[np.where(idx_tell)]
        return self.wave, self.flux


if __name__ == "__main__":
    filename = "/HD170740/RED_860/HD170740_w860_redl_20140915_O12.fits"
    sp = EdiblesSpectrum(filename)
    print("Barycentric Velocity is", sp.v_bary)
    print(sp.target)

    data = sp.getSpectrum(xmin=7661.5, xmax=7669)

    plt.plot(data[0], data[1], label="Geocentric")

    bary_data = sp.getSpectrum(xmin=7661.5, xmax=7669, bary=True)

    plt.plot(bary_data[0], bary_data[1], label="Barycentric")
    axes = plt.gca()
    # axes.set_xlim([7660, 7690])
    # axes.set_ylim([0, 2500])
    # plt.vlines((7661.5, 7669), 0, 2200, linestyles='dashed', colors='r')
    # plt.legend()
    # plt.show()

    filename = "/HD170740/RED_860/HD170740_w860_redl_20160613_O12.fits"
    sp = EdiblesSpectrum(filename)
    print("Barycentric Velocity is", sp.v_bary)
    print(sp.target)

    data = sp.getSpectrum(xmin=7661.5, xmax=7669)

    plt.plot(data[0], data[1], label="Geocentric")

    bary_data = sp.getSpectrum(xmin=7661.5, xmax=7669, bary=True)

    plt.plot(bary_data[0], bary_data[1], label="Barycentric")
    axes = plt.gca()
    # axes.set_xlim([7660, 7690])
    # axes.set_ylim([0, 2500])
    plt.vlines((7661.5, 7669), 0, 2200, linestyles="dashed", colors="r")
    plt.legend()
    plt.show()
