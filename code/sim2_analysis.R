#-----------------------------------------------------------#
#             Simulation analysis 2: virtual species        #
#-----------------------------------------------------------#

library("parallel")
library("doParallel")
library("pbapply")

# Load utils, functions, and define number of iterations
source("code/sim2_functions.R")
nsim <- 100
pboptions(type = "timer")

# Read data
spoly <- st_read(dsn="data/species_vdata.gpkg", layer="sampling_polygon")
wclim <- stack("data/species_stack.grd")
wgrid <- st_read(dsn="data/species_vdata.gpkg", layer="landscape_grid")

# Prepare parallelization
print(paste0("Process started with ", detectCores(), " cores."))
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Launch simulation
set.seed(1234)
sims <- pbreplicate(nsim, sim_species(wgrid, wclim, spoly), simplify=FALSE)

# We're done
stopCluster(cl)
rm("cl")
# write_csv(do.call(rbind, sims), "results/sim2_species/results_iberia.csv")
