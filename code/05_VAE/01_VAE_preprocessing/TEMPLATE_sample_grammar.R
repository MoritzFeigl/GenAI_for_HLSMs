# GenAI_para-ET-mHM project
# Sample Grammar for GenAI_para TEMPLATE FOR VSC

main_path <- "/gpfs/data/fs71468/GenAI_para_runs/"
setwd(main_path)
code_dir <- paste0(main_path, "code")
code_files <- list.files(paste0(code_dir, "/02_utility_functions/09_CFG_implementation"), full.names = TRUE)
for(file in code_files) source(file)
# Define variables and operators of the CFG
operators <- c('+', '-', '*', '/')
# Define grammar
grammar <- create_grammar(eq = "<eq> <op> <eq>,
                          <eq> <op> numeric,
                          <eq> <op> var,
                          <eq> <op> (<eq>),
                          var,
                          f(var),
                          f(<eq>),
                          (<eq>)^(<pm>numeric),
                          numeric",
                          op = paste(operators, collapse = ","),
                          pm = "+, -")
# Sample grammar
funs <- grammar_sampler(n = 36000000, grammar = grammar, max_depth = 8,
                        no_cores = 250, seed = 1304, save = TRUE,
                        file_name = "36mil_sampled_grammar")
