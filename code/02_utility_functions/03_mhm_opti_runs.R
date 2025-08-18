# GenAI_para-ET-mHM project
# experiment runs

# define run name to create environment automatically
experiment <- 1
run_name <- "exp1-run1"
main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"

# Environment setup ----------------------------------------------------------------------
# create run and result folders
run_dir <- paste0(main_path, "runs/", run_name)
result_dir <- paste0(main_path, "results/", run_name)
code_dir <- paste0(main_path, "code/")
data_dir <- paste0(main_path, "data/")
dir.create(run_dir, showWarnings = FALSE)
dir.create(result_dir, showWarnings = FALSE)
setwd(run_dir)

# load utils
source(paste0(code_dir, "02_utility_functions/01_GenAI_para_utility_functions.R"))
source(paste0(code_dir, "02_utility_functions/02_vae_generators.R"))

# define experiment configuration
source(paste0(code_dir, 
              "06_run_experiments/experiment_configs/experiment_", experiment, "_setup.R"))

# parallel setups is bound to number of complexes
n_parallel_setups <- ncomplexes

# create run environment
create_run_environment(run_dir, data_dir, code_dir, n_parallel_setups)

# define mhm output --> always aET for result analysis
state_and_fluxes[["aET"]] <- TRUE
for(setup in 1:n_parallel_setups){
  change_mhm_outputs(basins_dir = paste0(run_dir, "/Training_", setup), 
                     state_and_fluxes, 
                     time_step = "monthly")
}

# objective function ---------------------------------------------------------------------
objective_function <- objective_function_definer(run_dir, 
                                                 data_dir,
                                                 only_numerics = (experiment == 0),
                                                 parameter_nml = parameter_nml,
                                                 ind_of_num_paras,
                                                 loss_mode, 
                                                 latent_dim)

# previous optimizer state ---------------------------------------------------------------
# get previously computed optimizer states if available or else start initial compilation
if(!is.null(result_dir) & file.exists(paste0(result_dir, "/optimization_state.rds"))){
  opt_state <- readRDS(paste0(result_dir, "/optimization_state.rds"))
  obj <- readRDS(paste0(result_dir, "/last_sce_obj.rds"))
} else {
  counter <- 0
  obj <- list("iterations" = 0)
}

# initial population ---------------------------------------------------------------------
POPULATION <- sample_initial_popoulation(run_dir, experiment, obj,ncomplexes, NDIM, 
                                         vae_list, latent_dim, standard_parameters)

# # SCE optimization ---------------------------------------------------------------------
# Set SCE parameter
GenAI_para_dims <- ifelse(experiment %in% 1:2, latent_dim*tf_parameters, NA)
sceDefaults <- function(){
  list(ncomplex = ncomplexes,     ## number of complexes
       cce.iter = 2 * NDIM + 1,  ## number of iteration in inner loop (CCE algorithm) recommended number of CCE steps in Duan et al 1994:
       nPOINTS_COMPLEX = 2 * NDIM + 1, # points per complex
       fnscale = -1,            ## function scaling factor (set to -1 for maximisation)
       elitism = 1,            ## controls amount of weighting in sampling towards the better parameter sets
       initsample = "random",   ## sampling scheme for initial values -- "latin" or "random", only used if ini_pop = NULL
       reltol = 1e-5,          ## convergence threshold: relative improvement factor required in an SCE iteration
       tolsteps = 50,           ## number of iterations within reltol to confirm convergence
       maxit = numIter,          ## maximum number of iterations
       maxeval = Inf,          ## maximum number of function evaluations
       maxtime = Inf,          ## maximum duration of optimization in seconds
       returnpop = FALSE,      ## whether to return populations from all iterations
       trace = 0,              ## level of user feedback
       GenAI_para_dims = GenAI_para_dims,       # number of GenAI_para TF dimensions
       REPORT = 1,## number of iterations between reports when trace >= 1
       time_limit = 55) ## time limit for cce iterations in hours
}

# start running vae generator python scripts in background
if(experiment %in% 1:2) start_GenAI_para_vae(run_dir, code_dir, vae_list, ncomplexes)
system(paste0("chmod -R 777 ", run_dir))
system(paste0("chmod -R 777 ", result_dir))

if(obj$iterations < numIter){
  cat("Start/continue SCE optimization\n")
  obj <- SCEoptim(FUN = objective_function, 
                  par = start_point,
                  lower = xBounds$lower,
                  upper = xBounds$upper,
                  state_dir = result_dir,
                  restart_after_n_iter = restart_after_n_iter,
                  ini_pop = POPULATION)
}
if(obj$iterations < numIter) {
  system("sbatch 01_experiment_start.sh")
}

