import openmc
import numpy as np
from .mvng_source import mvng_source_diamonds
# needed to download cross sections on the fly
import openmc_data_downloader as odd
import matplotlib.pyplot as plt

from src import salt_volume, clif_density

salt_v = salt_volume

plt.rc('font', size=16) 

#   MATERIALS

# PbLi - natural - pure
pbli = openmc.Material(name="pbli")
pbli.add_element("Pb", 84.2, "ao")
pbli.add_element("Li", 15.2, "ao")
pbli.set_density("g/cm3", 11)

flibe = openmc.Material(name="flibe")
flibe.add_element("Li", 2.0, "ao")
flibe.add_element("Be", 1.0, "ao")
flibe.add_element("F", 4.0, "ao")
flibe.set_density("g/cm3", 1.94)

# lif-licl - natural - pure
clif = openmc.Material(name="clif")
clif.add_element("F", 0.5 * 0.305, "ao")
clif.add_element("Li", 0.5 * 0.305 + 0.5 * 0.695, "ao")
clif.add_element("Cl", 0.5 * 0.695, "ao")
clif.set_density("g/cm3", clif_density)

# FLiNaK - natural - pure
flinak = openmc.Material(name="flinak")
flinak.add_element("F", 50, "ao")
flinak.add_element("Li", 23.25, "ao")
flinak.add_element("Na", 5.75, "ao")
flinak.add_element("K", 21, "ao")
flinak.set_density("g/cm3", 2.020)

# thermal insulator
insulator = openmc.Material(name="insulator")
insulator.add_element("Al", 0.250, "ao")
insulator.add_element("Si", 0.125, "ao")
insulator.add_element("O", 0.625, "ao")
insulator.set_density("g/cm3", 0.128)

ss316 = openmc.Material(name="ss316")
ss316.add_element("Cr", 0.18)
ss316.add_element("Fe", 0.63)
ss316.add_element("Mn", 0.02)
ss316.add_element("Mo", 0.03)
ss316.add_element("Ni", 0.14)
ss316.set_density("g/cm3", 7.87)

inconel625 = openmc.Material(name="inconel625")
inconel625.add_element("C", 0.004852, "ao")
inconel625.add_element("Al", 0.008639, "ao")
inconel625.add_element("Si", 0.010374, "ao")
inconel625.add_element("P", 0.000281, "ao")
inconel625.add_element("S", 0.000272, "ao")
inconel625.add_element("Ti", 0.004869, "ao")
inconel625.add_element("Cr", 0.0243380, "ao")
inconel625.add_element("Mn", 0.005303, "ao")
inconel625.add_element("Fe", 0.052167, "ao")
inconel625.add_element("Co", 0.009887, "ao")
inconel625.add_element("Ni", 0.581642, "ao")
inconel625.add_element("Nb", 0.023124, "ao")
inconel625.add_element("Mo", 0.055210, "ao")
inconel625.set_density("g/cm3", 8.44)

air = openmc.Material(name="air")
air.add_element("C", 0.000150, "ao")
air.add_element("N", 0.784431, "ao")
air.add_element("O", 0.210748, "ao")
air.add_element("Ar", 0.004671, "ao")
air.set_density("g/cm3", 0.001205)

sparge = openmc.Material(name="sparge")
sparge.add_element("H", 0.03, "wo")
sparge.add_element("He", 0.97, "wo")
sparge.set_density("g/cm3", 0.0001589)

name_to_pretty = {
    clif.name: "ClLiF",
    flinak.name: "FLiNaK",
    flibe.name: "FLiBe",
    pbli.name: "PbLi",
    ss316.name: "SS316",
    inconel625.name: "Inconel625",
    air.name: "Air",
    sparge.name: "Sparge gas",
    insulator.name: "Insulator",
}

def make_model(breeder_material: openmc.Material, batches: int, particles: int):
    """Returns an openmc model for the BABY experiment with one type of breeder material

    Args:
        breeder_material (openmc.Material): the breeder material
        batches (int): number of batches
        particles (int): number of particles per batch

    Returns:
        openmc.Model: the openmc model for the BABY experiment
    """
    materials = openmc.Materials(
        [pbli, flibe, ss316, inconel625, air, sparge, clif, flinak, insulator]
    )

    materials.download_cross_section_data(
        libraries=["FENDL-3.1d"],
        set_OPENMC_CROSS_SECTIONS=True,
        particles=["neutron"],
    )

    # GEOMETRIC PARAMETERS --------------------------------------------------------
    # Numbered radii with low numbers correspond to inner radii and high numbers
    # move concentrically outward. Numbered "z"-values go from lowest to highest
    # both in number and position.

    # neutron generator
    gen_t = 0.25  # generator wall thickness
    gen_ri = 4.75  # generator inner radius
    gen_xo_min = -36.5  # generator legth behind the Zr-T target
    gen_xo_max = 13.5  # generator length in front of the Zr-T target

    cru_t = 0.15  # crucible thickness
    cru_socket_ri = 0.238  # crucible socket for heater inner radius
    cru_socket_h = 9.17  # socket inner height
    cru_ri = 2.1  # crucible inner radius
    cru_h = 10.17  # crucible height

    # crucible position
    z_bottom_crucible = -3.0

    # salt_mass = 190  # g salt mass measured by Weiyue
    # salt_v = salt_mass / clif.density  # use the volume of ClLiF for all breeders
    salt_h = (salt_v / np.pi - (cru_h - cru_socket_h) * cru_ri**2) / (
        cru_ri**2 - (cru_socket_ri + cru_t) ** 2
    )
    print(f"Breeder volume: {salt_v:.2e} cm3")
    salt_h += cru_h - cru_socket_h  # salt height in the crucible

    # GEOMETY

    # surfaces

    # neutron generator
    cx_1 = openmc.XCylinder(y0=0, z0=0, r=gen_ri)  # generator inner body
    # generator body + thickness
    cx_2 = openmc.XCylinder(y0=0, z0=0, r=gen_ri + gen_t)
    px_3 = openmc.XPlane(x0=gen_xo_min)  # generator bottom face (far left)
    # generator bottom face + wall thickness
    px_4 = openmc.XPlane(x0=gen_xo_min + gen_t)
    px_5 = openmc.XPlane(x0=gen_xo_max - gen_t)  # generator top face
    # generator top face + thickness (far right)
    px_6 = openmc.XPlane(x0=gen_xo_max)

    # second generator
    cx_30 = openmc.XCylinder(y0=24.1, z0=0, r=gen_ri)  # generator inner body
    # generator body + thickness
    cx_31 = openmc.XCylinder(y0=24.1, z0=0, r=gen_ri + gen_t)

    # crucible and vessel

    y0 = 12.7  # 10.5

    # crucible inner socket
    cz_7 = openmc.ZCylinder(x0=0, y0=y0, r=cru_socket_ri)
    # crucible inner socket + thickness
    cz_8 = openmc.ZCylinder(x0=0, y0=y0, r=cru_socket_ri + cru_t)
    cz_9 = openmc.ZCylinder(x0=0, y0=y0, r=cru_ri)  # crucible inner body
    # crucible inner body + thickness
    cz_10 = openmc.ZCylinder(x0=0, y0=y0, r=cru_ri + cru_t)
    pz_14 = openmc.ZPlane(z0=z_bottom_crucible - 2.0)  # bottom insulator
    pz_15 = openmc.ZPlane(z0=z_bottom_crucible)  # bottom crucible
    # bottom crucible + thickness
    pz_16 = openmc.ZPlane(z0=z_bottom_crucible + cru_t)
    pz_17 = openmc.ZPlane(
        z0=z_bottom_crucible + cru_t + cru_h - cru_socket_h - cru_t
    )  # crucible socket bottom
    # crucible socket bottom + thickness
    pz_18 = openmc.ZPlane(z0=z_bottom_crucible + cru_t + cru_h - cru_socket_h)
    pz_19 = openmc.ZPlane(z0=z_bottom_crucible + cru_t +
                          salt_h)  # flibe free surface
    pz_20 = openmc.ZPlane(z0=z_bottom_crucible + cru_t + cru_h)  # top cucible
    pz_21 = openmc.ZPlane(
        z0=z_bottom_crucible + cru_t + cru_h + cru_t
    )  # top crucible + thickness
    pz_22 = openmc.ZPlane(
        z0=z_bottom_crucible + cru_t + cru_h + cru_t + 2
    )  # top insulator

    # insulator
    cz_24 = openmc.ZCylinder(x0=0, y0=y0, r=cru_ri + cru_t + 2.5)
    cz_25 = openmc.ZCylinder(x0=0, y0=y0, r=cru_ri + cru_t + 2.8)

    # outer region
    so_999 = openmc.Sphere(x0=0, y0=0, z0=0, r=500.0, boundary_type="vacuum")

    # regions
    # neutron generator 1
    region_1 = +px_4 & -px_5 & -cx_1  # inner void
    region_2 = +px_3 & -px_6 & -cx_2 & (-px_4 | +px_5 | +cx_1)  # ss shell
    # neutron generator 2
    region_30 = +px_4 & -px_5 & -cx_30  # inner void
    region_31 = +px_3 & -px_6 & -cx_31 & (-px_4 | +px_5 | +cx_30)  # ss shell

    # crucible, flibe and vessel
    region_3 = +pz_16 & -pz_19 & -cz_9 & (-pz_17 | +cz_8)  # salt
    region_4 = +pz_19 & -pz_20 & -cz_9 & +cz_8  # sparge gas in crucible
    region_5 = (
        +pz_15
        & -pz_21
        & -cz_10
        & (-pz_16 | +pz_20 | +cz_9 | (+pz_17 & -cz_8))
        & (-pz_18 | +cz_7)
    )

    region_8 = +pz_14 & -pz_22 & +cz_10 & -cz_24
    region_9 = +pz_14 & -pz_22 & +cz_24 & -cz_25

    region_999 = -so_999 & ~(
        region_1
        | region_2
        | region_30
        | region_31
        | region_3
        | region_4
        | region_5
        | region_8
        | region_9
    )

    # cells

    # neutron generator
    cell_1 = openmc.Cell(cell_id=1, fill=None,
                         region=region_1)  # generator void
    cell_2 = openmc.Cell(cell_id=2, fill=ss316,
                         region=region_2)  # generator shell
    cell_30 = openmc.Cell(cell_id=30, fill=None,
                          region=region_30)  # generator void
    cell_31 = openmc.Cell(cell_id=31, fill=ss316,
                          region=region_31)  # generator shell
    # salt
    cell_3 = openmc.Cell(cell_id=3, fill=breeder_material,
                         region=region_3)  # salt
    # sparge gas in crucible
    cell_4 = openmc.Cell(cell_id=4, fill=sparge, region=region_4)
    cell_5 = openmc.Cell(cell_id=5, fill=inconel625,
                         region=region_5)  # crucible shell

    cell_8 = openmc.Cell(cell_id=8, fill=insulator,
                         region=region_8)  # insulator
    cell_9 = openmc.Cell(cell_id=9, fill=insulator,
                         region=region_9)  # insulator

    cell_999 = openmc.Cell(cell_id=999, fill=air,
                           region=region_999)  # outer sphere

    universe = openmc.Universe(
        name="universe",
        cells=[
            cell_1,
            cell_2,
            cell_30,
            cell_31,
            cell_3,
            cell_4,
            cell_5,
            cell_8,
            cell_9,
            cell_999,
        ],
    )
    # plot geometry with materials
    mat_to_colour = {
        breeder_material: (30, 12, 245),
        ss316: (74, 73, 108),
        inconel625: (181, 38, 24),
        air: (190, 253, 254),
        sparge: (188, 252, 143),
        insulator: (254, 255, 145),
    }
    universe.plot(
        origin=(0,12.7,2),
        width=(37, 16),
        pixels=int(7e5),
        basis='yz',
        color_by='material',
        colors=mat_to_colour,
        outline=True,
        legend=True,
        legend_kwargs={"fontsize": 16, "framealpha": 1}
        )
    # replace labels in legend
    for text in plt.gca().get_legend().get_texts():
        pretty_label = name_to_pretty[text.get_text()]
        text.set_text(pretty_label)

    plt.tight_layout()
    for ext in ["png", "svg", "pdf"]:
        plt.savefig(f"geometry_{breeder_material.name}.{ext}")
    geometry = openmc.Geometry(universe)
    geometry.merge_surfaces = True

    # TALLIES

    # cell filters
    breeder_filter = openmc.CellFilter([cell_3])

    # energy filters

    # tallies

    # cell tally - salt tbr
    tally1 = openmc.Tally(1, "salt_cell_tbr")
    tally1.filters = [breeder_filter]
    tally1.scores = ["(n,Xt)"]
    tally1.nuclides = ["Li6", "Li7", "F19"]
    #
    tally2 = openmc.Tally(2, "salt_cell_tbr2")
    tally2.filters = [breeder_filter]
    tally2.scores = ["(n,Xt)"]

    tallies = openmc.Tallies([tally1, tally2])

    # SETTINGS

    # source definition

    # mvng source characterized via diamond detectors
    source1 = mvng_source_diamonds(center=(0, 0, 0), reference_uvw=(0, 1, 0))
    source2 = mvng_source_diamonds(
        center=(0, 24.1, 0), reference_uvw=(0, 1, 0))

    my_source = []
    for s in source1:
        s.strength *= 0.233
        my_source.append(s)
    for s in source2:
        s.strength *= 0.767
        my_source.append(s)

    # settings
    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.source = my_source
    settings.batches = batches
    settings.particles = int(particles)
    settings.output = {"tallies": False}

    model = openmc.Model(
        materials=materials, geometry=geometry, settings=settings, tallies=tallies
    )
    return model

def main(batches: int = 100, particles: int = int(1e7)):
    """Main function running the openmc model for four different breeders

    Args:
        batches (int, optional): Number of batches. Defaults to 100.
        particles (int, optional): Number of particles per batch. Defaults to int(1e7).
    """
    for breeder_material in [pbli, flibe, clif, flinak]:
        model = make_model(breeder_material, batches=batches,
                           particles=particles)
        model.run(
            threads=16,
            cwd=breeder_material.name,
        )


if __name__ == "__main__":
    main(batches=100, particles=int(1e6))
