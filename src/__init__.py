from pathlib import Path
import inspect

def absolute_path(filename: str, level=1):
    """Returns the absolute path of a file. Based on a relative path.

    Args:
        filename (str): the relative path to the file
        level (int, optional): Level in the file call. 0 corresponds to
            the file where absolute_path is defined, 1 correspond to the
            file calling this function, 2 corresponds to the parent of
            the file calling this function. Defaults to 1.

    Returns:
        str: the absolute path of the file
    """
    caller_frame = inspect.stack()[level]
    return str(Path(caller_frame.filename).parent) + "/" + filename

def get_cllif_density(temp, LiCl_frac=0.695, cl37_enr=0.2424):
    """Calculates ClLiF density (g/cc) given a temperature (Celsius),
    LiCl molar fraction, and Cl-37 enrichment by calculating the molar volume of
    LiCl and LiF from published correlations, 
    using those values to calculate the mixture molar volume, and dividing the
    ClLiF molecular weight by the ClLiF molar volume to get density.

    Args:
        temp (float): temperature in Celsius
        LiCl_frac (float, optional): LiCl molar fraction. Defaults to 0.695.
        cl37_enr (float, optional): Cl-37 enrichment. Defaults to 0.2424.

    Returns:
        float: density of ClLiF in g/cc
    """

    M_Cl35 = 34.969 #amu
    M_Cl37 = 36.9659 #amu
    M_Li = 6.941 #amu

    # Calculate LiCl density
    # The temperature density correlation for LiCl works up to 781 C
    # Source: http://moltensalt.org/references/static/downloads/pdf/element-salt-densities.pdf
    rho_m_LiCl = 1.502 #g/cc
    k_LiCl = 0.000432
    t_m_LiCl = 610
    rho_LiCl = rho_m_LiCl - k_LiCl*(temp - t_m_LiCl)

    # Calculate LiF density
    # This temperature density correlation for LiF works up to 1047 C
    # Source: http://moltensalt.org/references/static/downloads/pdf/element-salt-densities.pdf
    rho_m_LiF = 1.81 #g/cc
    k_LiF = 0.000490
    t_m_LiF = 848.2
    rho_LiF = rho_m_LiF - k_LiF*(temp - t_m_LiF)

    ## Calculate molar volume
    # Molar Mass of LiCl
    M_LiCl = M_Cl35 * (1 - cl37_enr) + M_Cl37 * cl37_enr + M_Li #g/mol
    V_LiCl = M_LiCl/rho_LiCl

    M_LiF = 25.939 #g/mol
    V_LiF = M_LiF/rho_LiF

    # molar volume of mixture
    V_m = LiCl_frac*V_LiCl + (1-LiCl_frac)*V_LiF
    # molar mass of mixture
    M_m = LiCl_frac*M_LiCl + (1-LiCl_frac)*M_LiF

    density_cllif = M_m/V_m
    return density_cllif

salt_mass = 190  # g
salt_temperature = 700 # C
clif_density = get_cllif_density(salt_temperature, LiCl_frac=0.695, cl37_enr=0.2424)
salt_volume = salt_mass / clif_density  # use the volume of ClLiF for all breeders

from .openmc_model.openmc_model import main

