# building the MIT-VaultLab neutron generator
# angular and energy distribution

import numpy as np
import h5py
import openmc

from pathlib import Path
from src import absolute_path
path = absolute_path('mvng_source_diamonds.h5')

mvng_source = h5py.File(path, 'r')

angles = ["0", "15", "30", "45", "60", "75", "90", "105", "120", "135", "150"]


def mvng_source_diamonds(center=(0, 0, 0), reference_uvw=(0, 0, 1)):
    '''method for building the MIT-VaultLab neutron generator in OpenMC
    with data tabulated from John Ball and Shon Mackie characterization
    via diamond detectors

    Parameters
    ----------
    center : coordinate position of the source (it is a point source)

    reference_uvw : direction for the polar angle (tuple or list of versors)
    it is the same for the openmc.PolarAzimuthal class
    more specifically, polar angle = 0 is the direction of the D accelerator
    towards the Zr-T target
    '''

    # energy values
    energies = mvng_source['values/table']['Energy (MeV)'] * 1e6
    # angular bins in [0, pi)
    pbins = np.cos([np.deg2rad(float(a)) for a in angles] + [np.pi])
    spectra = [mvng_source['values/table'][col] for col in angles]

    # yield values for strengths
    yields = np.sum(spectra, axis=-1) * np.diff(pbins)
    yields /= np.sum(yields)

    # azimuthal values
    phi = openmc.stats.Uniform(a=0, b=2*np.pi)

    all_sources = []
    for i, angle in enumerate(pbins[:-1]):

        mu = openmc.stats.Uniform(a=pbins[i+1], b=pbins[i])

        space = openmc.stats.Point(center)
        angle = openmc.stats.PolarAzimuthal(
            mu=mu, phi=phi, reference_uvw=reference_uvw)
        energy = openmc.stats.Tabular(
            energies, spectra[i], interpolation='linear-linear')
        strength = yields[i]

        my_source = openmc.Source(
            space=space, angle=angle, energy=energy, strength=strength, particle='neutron')

        all_sources.append(my_source)

    return all_sources
