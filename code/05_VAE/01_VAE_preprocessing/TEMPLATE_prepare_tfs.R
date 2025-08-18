# Prepare transfer functions for GenAI_para
# GenAI_para-mHM project
# Moritz Feigl, Aug 2020


main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
code_dir <- paste0(main_path, "code")
code_files <- list.files(paste0(code_dir, "/02_utility_functions/09_CFG_implementation"), full.names = TRUE)
for(file in code_files) source(file)
setwd(paste0(main_path, "data/07_VAE_data/03_prepared_TFs"))
# Sample variables and numerics
# use common functions that map from R -> R
f <- c('sqrt', 'exp', 'log10', 'log', 
       'sin', 'sinh', 
       'cos', 'cosh', 
       'tan', 'tanh', 'atan',
       'abs')
base <- c("dem", "aspect_sin", "aspect_cos", "slope", "bd", "sand", "clay")
base_plus <- c("lai", "map", "mat", "mat_range")
mhm <- c("ThetaS", "KSat", "vGenu_n")
lai <- "lai"

vae_1 <- c(base)
vae_2 <- c(base, base_plus)
vae_3 <- c(base, mhm)
vae_4 <- c(base, base_plus, mhm)
vae_5 <- c(lai)
# specific selection
vae_version <- 1
var <- get(paste0("vae_", vae_version))
numbers <- seq(0.01, 100, 0.01)
variable_list <- list("var" = var,
                      "numeric" = numbers,
                      "f" = f)

variable_input(functions = "dummy_file",
               variable_list = variable_list,
               single_var_to_remove = "numeric",
               necessary_var_to_be_included = "var",
               n_iter = 20, no_cores = 110, seed = 1304,
               file_name = "dummy_functions")

simplify_functions(files = "dummy_files",
                   no_cores = 110)