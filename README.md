# Spatial simulation and models R-INLA

Simulation of three different sampling scenarios: (1) Random sampling (2) Preferential sampling and (3) Lattice data sampling. For each of them, biomass/catches (response variable), bathymetry and effort (explanatory variables) were simulated.

The repository contains the development of three different simulations. All of them are pretended to simulate the catches per unit effort (CPUE) in a spatial area during a period of 1 year. The three cases are the following:

1. Random sampling: simulation of an oceanographic survey with georeferenced data (independence in sampling). In this firts case, the response variable has more than 30 % zeros. 

2. Preferential sampling: simulation of fisheries with georeferenced data (loss of independence in sampling).

3. Areal data: simulation of the preferential sampling with spatial lattice data.

The following variables are simulated for each scenario:

(a) Biomass/catches: response variable measured in tons. It is a quantitative, continuous and positive variable. 

b) Effort: quantitative and continuous explanatory variable, measuring the time during which the gear is active.

c) Bathymetry: quantitative and continuous explanatory variable, measures depth in meters.

Once the simulation is completed, model fitting is performed in R-INLA for the respective scenarios. 
