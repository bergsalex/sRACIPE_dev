rm(list = ls())
library(reshape2)
library(ggplot2)
library(Rtsne)
library(gplots)
library(MASS)
library(RColorBrewer)
library(htmlwidgets)
library(d3heatmap)
library(sRACIPEv03)
working_directory <- getwd()
config_file = "inputs/sRACIPE.cfg"
configuration <- sRACIPE_load_configuration(config_file = config_file)

topology_file <- "/Users/koharv/Documents/GitHub/sRACIPE_dev_JAX/inputs/test2.tpo"
topology <- sRACIPE_load_topology(topology_file = topology_file)

results_directory <- ifelse(!dir.exists(file.path(working_directory, "results")), dir.create(file.path(working_directory, "results")), file.path(working_directory, "results"))

output_file <- sRACIPE_simulate_GRN(topology_file = topology_file, PARAMETERS_FILE = 1, MAX_NOISE = 0, NOISE_LEVELS = 1, NUM_MODELS = 10)

plot_data_stochastic(output_file, plot_filename = topology$filename)

#plot_network(topology)

#output_file <- sRACIPE_RK_deterministic(topology_file = topology_file)

#output_file <- sRACIPE_RK_adaptive_deterministic(topology_file = topology_file)

#output_file <- sRACIPE_stochastic(topology_file = topology_file)

output_file <- sRACIPE_RK_deterministic_knockout(topology_file = topology_file)

output_file <- sRACIPE_RK_deterministic_knockout(topology_file = topology_file, KNOCKOUT = "GCM1",NUM_MODELS = 2000)


plot_data_deterministic(output_file, plot_filename = topology$filename)

#plot_data_stochastic(output_file, plot_filename=topology$filename, topology_df=topology, config = configuration, bin_count=40)

plot_interactive_heatmap(output_file, plot_filename = topology$filename)
