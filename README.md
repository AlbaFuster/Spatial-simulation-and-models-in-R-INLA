# Spatial simulation and random sampling models R-INLA

Simulation of three different spatial scenarios: (1) Random sampling (2) Preferential sampling and (3) Areal data sampling. For each of them, abundance/catches (response variable), bathymetry and effort (explanatory variables) were simulated.

The repository contains the development of three different simulations. All of them are pretended to simulate the catches per unit effort in a spatial area during a period of 1 year. The three cases are the following:

1. Random sampling: simulation of an oceanographic campaign with georeferenced data (independence in sampling). Two cases are included: (1) Simulation with less than 30 % zeros (1.1) Simulation with more than 30 % zeros. 

2. Preferential sampling: simulation of fisheries with georeferenced data (loss of independence in sampling).

3. Areal data: simulation of the preferential sampling with spatial lattice data.

The following variables are simulated for each scenario:

(a) Biomass/catches: response variable measured in tons.

b) Effort: quantitative and continuous explanatory variable, measuring the time during which the gear is active.

c) Bathymetry: quantitative and continuous explanatory variable, measures depth in meters.

Once the simulation was carried out, the models for scenario (1) random sampling has been adjusted. In this case, the following two blocks of models have been adjusted: 

I.  Models without zeros: the response variables catches and CPUE have been transformed by adding 0.1

II. Models with zeros or hurdle models: zeros are modeled through a Bernouilli process. 
