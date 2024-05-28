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
library("FNN")
library("gstat")
library("CAST")
library("RandomForestsGLS")


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
sim10_partial <- lapply(1:nsim, sim_fields, range=10, scenario="partial")
sim10_partial <- do.call(rbind, sim10_partial)
write.csv(sim10_partial, paste0(pathroot, "results/simulations/sim10_partial.csv"), row.names=F)

# range=10 - complete
print(Sys.time())
sim10_complete <- lapply(1:nsim, sim_fields, range=10, scenario="complete")
sim10_complete <- do.call(rbind, sim10_complete)
write.csv(sim10_complete, paste0(pathroot, "results/simulations/sim10_complete.csv"), row.names=F)

# range=10 - proxyonly
print(Sys.time())
sim10_proxyonly <- lapply(1:nsim, sim_fields, range=10, scenario="proxyonly")
sim10_proxyonly <- do.call(rbind, sim10_proxyonly)
write.csv(sim10_proxyonly, paste0(pathroot, "results/simulations/sim10_proxyonly.csv"), row.names=F)

# range=10 - autocor
print(Sys.time())
sim10_autocor <- lapply(1:nsim, sim_fields, range=10, scenario="autocor")
sim10_autocor <- do.call(rbind, sim10_autocor)
write.csv(sim10_autocor, paste0(pathroot, "results/simulations/sim10_autocor.csv"), row.names=F)





# range=40 - partial
print(Sys.time())
sim40_partial <- lapply(1:nsim, sim_fields, range=40, scenario="partial")
sim40_partial <- do.call(rbind, sim40_partial)
write.csv(sim40_partial, paste0(pathroot, "results/simulations/sim40_partial.csv"), row.names=F)

# range=40 - complete
print(Sys.time())
sim40_complete <- lapply(1:nsim, sim_fields, range=40, scenario="complete")
sim40_complete <- do.call(rbind, sim40_complete)
write.csv(sim40_complete, paste0(pathroot, "results/simulations/sim40_complete.csv"), row.names=F)

# range=40 - proxyonly
print(Sys.time())
sim40_proxyonly <- lapply(1:nsim, sim_fields, range=40, scenario="proxyonly")
sim40_proxyonly <- do.call(rbind, sim40_proxyonly)
write.csv(sim40_proxyonly, paste0(pathroot, "results/simulations/sim40_proxyonly.csv"), row.names=F)

# range=40 - autocor
print(Sys.time())
sim40_autocor <- lapply(1:nsim, sim_fields, range=40, scenario="autocor")
sim40_autocor <- do.call(rbind, sim40_autocor)
write.csv(sim40_autocor, paste0(pathroot, "results/simulations/sim40_autocor.csv"), row.names=F)


# We're done
stopCluster(cl)
rm("cl")
