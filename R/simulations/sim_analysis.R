#------------------------------------------------------------------------------#
#                           Simulation analysis in HPC                         #
#------------------------------------------------------------------------------#

# HPC
pathroot <- "****"

# Packages
library("parallel")
library("doParallel")
library("pbapply")
library("caret")
library("randomForest")
library("sf")
library("raster")
library("dplyr")
library("tidyr")

# Requirements local packages
library("FNN")

# Local packages
# library("gstat")
library("gstat", lib.loc = paste0(pathroot, "packages/"))
# library("CAST")
library("CAST", lib.loc = paste0(pathroot, "packages/"))
# library("RandomForestsGLS")
library("RandomForestsGLS", lib.loc = paste0(pathroot, "packages/"))

# other dependencies
library("twosamples", lib.loc = paste0(pathroot, "packages/"))

# Load utils, functions, and define number of iterations
source(paste0(pathroot, "R/simulations/sim_functions.R"))
source(paste0(pathroot, "R/simulations/sim_utils.R"))
nsim <- 100
pboptions(type = "timer")

# Prepare parallelization
cl <- makeCluster(8)
registerDoParallel(cl)

# range=10 - partial
print(Sys.time())
set.seed(1234)
sim10_partial <- pbreplicate(nsim, sim_fields(range=10, simtype="partial"), simplify=FALSE)
sim10_partial <- do.call(rbind, sim10_partial)
write.csv(sim10_partial, paste0(pathroot, "results/simulations/sim10_partial.csv"), row.names=F)

# range=10 - complete
print(Sys.time())
set.seed(1234)
sim10_complete <- pbreplicate(nsim, sim_fields(range=10, simtype="complete"), simplify=FALSE)
sim10_complete <- do.call(rbind, sim10_complete)
write.csv(sim10_complete, paste0(pathroot, "results/simulations/sim10_complete.csv"), row.names=F)

# range=40 - partial
print(Sys.time())
set.seed(1234)
sim40_partial <- pbreplicate(nsim, sim_fields(range=40, simtype="partial"), simplify=FALSE)
sim40_partial <- do.call(rbind, sim40_partial)
write.csv(sim40_partial, paste0(pathroot, "results/simulations/sim40_partial.csv"), row.names=F)

# range=40 - complete
print(Sys.time())
set.seed(1234)
sim40_complete <- pbreplicate(nsim, sim_fields(range=40, simtype="complete"), simplify=FALSE)
sim40_complete <- do.call(rbind, sim40_complete)
write.csv(sim40_complete, paste0(pathroot, "results/simulations/sim40_complete.csv"), row.names=F)

# We're done
stopCluster(cl)
rm("cl")
