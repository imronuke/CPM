from math import log10

import numpy as np

import openmc
import openmc.mgxs

###############################################################################
# Create multigroup data

# Instantiate the energy group data
groups = openmc.mgxs.EnergyGroups(group_edges=[
    1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6])

# Instantiate the 7-group (C5G7) cross section data
uo2_xsdata = openmc.XSdata('UO2', groups)
uo2_xsdata.order = 0
uo2_xsdata.set_total(
    [1.7795E-01,3.2980E-01,4.8039E-01,5.5437E-01,3.1180E-01,3.9517E-01,5.6441E-01])
uo2_xsdata.set_absorption([8.0248E-03, 3.7174E-03, 2.6769E-02, 9.6236E-02,
                           3.0020E-02, 1.1126E-01, 2.8278E-01])
scatter_matrix = np.array(
    [[[1.2754E-01,4.2378E-02,9.4000E-06,0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00],
      [0.0000E+00,3.2446E-01,1.6314E-03,0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00],
      [0.0000E+00,0.0000E+00,4.5094E-01,2.6792E-03,0.0000E+00,0.0000E+00,0.0000E+00],
      [0.0000E+00,0.0000E+00,0.0000E+00,4.5257E-01,5.5664E-03,0.0000E+00,0.0000E+00],
      [0.0000E+00,0.0000E+00,0.0000E+00,1.2530E-04,2.7140E-01,1.0255E-02,0.0000E+00],
      [0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.2968E-03,2.6580E-01,1.6809E-02],
      [0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,8.5458E-03,2.7308E-01]]])
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
uo2_xsdata.set_scatter_matrix(scatter_matrix)
uo2_xsdata.set_fission([7.21206E-03, 8.19301E-04, 6.45320E-03,
                        1.85648E-02, 1.78084E-02, 8.30348E-02,
                        2.16004E-01])
uo2_xsdata.set_nu_fission([2.0060E-02,2.0273E-03,1.5706E-02,4.5183E-02,4.3342E-02,2.0209E-01,5.2571E-01
])
uo2_xsdata.set_chi([5.8791E-01, 4.1176E-01, 3.3906E-04, 1.1761E-07, 0.0000E+00,
                    0.0000E+00, 0.0000E+00])

h2o_xsdata = openmc.XSdata('LWTR', groups)
h2o_xsdata.order = 0
h2o_xsdata.set_total([1.5921E-01,4.1297E-01,5.9031E-01,5.8435E-01,7.1800E-01,1.2544E+00,2.6504E+00])

h2o_xsdata.set_absorption([6.0105E-04, 1.5793E-05, 3.3716E-04,
                           1.9406E-03, 5.7416E-03, 1.5001E-02,
                           3.7239E-02])
scatter_matrix = np.array(
    [[[4.4478E-02,1.1340E-01,7.2350E-04,3.7000E-06,1.0000E-07,0.0000E+00,0.0000E+00],
      [0.0000E+00,2.8233E-01,1.2994E-01,6.2340E-04,4.8000E-05,7.4000E-06,1.0000E-06],
      [0.0000E+00,0.0000E+00,3.4526E-01,2.2457E-01,1.6999E-02,2.6443E-03,5.0340E-04],
      [0.0000E+00,0.0000E+00,0.0000E+00,9.1028E-02,4.1551E-01,6.3732E-02,1.2139E-02],
      [0.0000E+00,0.0000E+00,0.0000E+00,7.1400E-05,1.3914E-01,5.1182E-01,6.1229E-02],
      [0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,2.2157E-03,6.9991E-01,5.3732E-01],
      [0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.3244E-01,2.4807E+00]]])
scatter_matrix = np.rollaxis(scatter_matrix, 0, 3)
h2o_xsdata.set_scatter_matrix(scatter_matrix)

mg_cross_sections_file = openmc.MGXSLibrary(groups)
mg_cross_sections_file.add_xsdatas([uo2_xsdata, h2o_xsdata])
mg_cross_sections_file.export_to_hdf5()

###############################################################################
# Create materials for the problem

# Instantiate some Macroscopic Data
uo2_data = openmc.Macroscopic('UO2')
h2o_data = openmc.Macroscopic('LWTR')

# Instantiate some Materials and register the appropriate Macroscopic objects
uo2 = openmc.Material(name='UO2 fuel')
uo2.set_density('macro', 1.0)
uo2.add_macroscopic(uo2_data)

water = openmc.Material(name='Water')
water.set_density('macro', 1.0)
water.add_macroscopic(h2o_data)

# Instantiate a Materials collection and export to XML
materials_file = openmc.Materials([uo2, water])
materials_file.cross_sections = "mgxs.h5"
materials_file.export_to_xml()

###############################################################################
# Define problem geometry

# Create a surface for the fuel outer radius
fuel_rad = openmc.ZCylinder(r=0.40, name='Fuel_rad')

# Create a region represented as the inside of a rectangular prism
water_rad = openmc.ZCylinder(r=0.70, name='water_rad', boundary_type='white')

# Instantiate Cells
fuel = openmc.Cell(fill=uo2, region=-fuel_rad, name='fuel')
moderator = openmc.Cell(fill=water, region=+fuel_rad & -water_rad, name='moderator')

# Create a geometry with the two cells and export to XML
geometry = openmc.Geometry([fuel, moderator])
geometry.export_to_xml()

###############################################################################
# Define problem settings

# Instantiate a Settings object, set all runtime parameters, and export to XML
settings = openmc.Settings()
settings.energy_mode = "multi-group"
settings.batches = 120
settings.inactive = 20
settings.particles = 75000

# Create an initial uniform spatial source distribution over fissionable zones
pitch = 1.4
lower_left = (-pitch/2, -pitch/2, -1)
upper_right = (pitch/2, pitch/2, 1)
uniform_dist = openmc.stats.Box(lower_left, upper_right)
settings.source = openmc.IndependentSource(
    space=uniform_dist, constraints={'fissionable': True})
settings.export_to_xml()

###############################################################################
#                   Exporting to OpenMC plots.xml file
###############################################################################

plot = openmc.Plot(plot_id=1)
plot.origin = [0, 0, 0]
plot.width = [4, 4]
plot.pixels = [400, 400]
plot.color_by = 'material'

# Instantiate a Plots collection and export to XML
plot_file = openmc.Plots([plot])
plot_file.export_to_xml()
