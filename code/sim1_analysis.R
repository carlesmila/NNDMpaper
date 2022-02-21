#-----------------------------------------------------------#
#            Simulation analyses 1: Random fields           #
#-----------------------------------------------------------#

library("parallel")
library("doParallel")
library("pbapply")

# Load utils, functions, and define number of iterations
source("code/sim1_functions.R")
nsim <- 100
pboptions(type = "timer")

# Prepare parallelization
print(paste0("Process started with ", detectCores(), " cores."))
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# range=1
print(Sys.time())
set.seed(1234)
sim1 <- pbreplicate(nsim, sim_fields(range=1), simplify=FALSE)
sim1 <- do.call(rbind, sim1)
# write.csv(sim1, paste0("results/sim1_fields/sim1.csv"), row.names=F)

# range=10
print(Sys.time())
set.seed(1234)
sim10 <- pbreplicate(nsim, sim_fields(range=10), simplify=FALSE)
sim10 <- do.call(rbind, sim10)
# write.csv(sim10, paste0("results/sim1_fields/sim10.csv"), row.names=F)

# range=20
print(Sys.time())
set.seed(1234)
sim20 <- pbreplicate(nsim, sim_fields(range=20), simplify=FALSE)
sim20 <- do.call(rbind, sim20)
# write.csv(sim20, paste0("results/sim1_fields/sim20.csv"), row.names=F)

# range=30
print(Sys.time())
set.seed(1234)
sim30 <- pbreplicate(nsim, sim_fields(range=30), simplify=FALSE)
sim30 <- do.call(rbind, sim30)
# write.csv(sim30, paste0("results/sim1_fields/sim30.csv"), row.names=F)

# range=40
print(Sys.time())
set.seed(1234)
sim40 <- pbreplicate(nsim, sim_fields(range=40), simplify=FALSE)
sim40 <- do.call(rbind, sim40)
# write.csv(sim40, paste0("results/sim1_fields/sim40.csv"), row.names=F)

# We're done
stopCluster(cl)
rm("cl")
