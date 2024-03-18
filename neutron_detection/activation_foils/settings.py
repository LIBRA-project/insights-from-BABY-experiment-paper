import pint
import numpy as np

ureg = pint.UnitRegistry()

geometric_efficiency = 0.5
nal_gamma_efficiency = 0.344917296922981
branching_ratio = 0.999


# source: https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WPJET3PR17_16948_submitted.pdf
# See Figure 3
Nb93_n_2n_Nb92m_cross_section_at_14Mev = 0.46 * ureg.barn

# source https://gammaray.inl.gov/SiteAssets/catalogs/ge/pdf/nb92m.pdf
Nb92m_half_life = 10.25 * ureg.day
Nb92m_decay_constant = np.log(2) / Nb92m_half_life
