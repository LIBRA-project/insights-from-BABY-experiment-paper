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