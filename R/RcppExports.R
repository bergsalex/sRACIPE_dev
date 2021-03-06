# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

interaction_reader <- function(gene_interaction, filepath, filename, number_gene) {
    .Call('_sRACIPEv03_interaction_reader', PACKAGE = 'sRACIPEv03', gene_interaction, filepath, filename, number_gene)
}

simulate_GRN <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, D_max, D_shot_scaling, GENE_NOISE_SCALING, file_writing_interval, D_levels, D_scaling, output_precision, ANNEALING, CONSTANT_NOISE, INITIAL_CONDITIONS, filename, parameters_file, readIC) {
    .Call('_sRACIPEv03_simulate_GRN', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, D_max, D_shot_scaling, GENE_NOISE_SCALING, file_writing_interval, D_levels, D_scaling, output_precision, ANNEALING, CONSTANT_NOISE, INITIAL_CONDITIONS, filename, parameters_file, readIC)
}

multiGeneCircuit_EM_uniform_Darray_annealing <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, D_max, D_shot_scaling, GENE_NOISE_SCALING, file_writing_interval, D_levels, D_scaling, output_precision, ANNEALING, CONSTANT_NOISE, INITIAL_CONDITIONS, filename) {
    .Call('_sRACIPEv03_multiGeneCircuit_EM_uniform_Darray_annealing', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, D_max, D_shot_scaling, GENE_NOISE_SCALING, file_writing_interval, D_levels, D_scaling, output_precision, ANNEALING, CONSTANT_NOISE, INITIAL_CONDITIONS, filename)
}

multiGeneCircuit_EM_uniform_Darray_multiprint <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, D_max, D_shot_scaling, GENE_NOISE_SCALING, file_writing_interval, D_levels, D_scaling, output_precision, ANNEALING, CONSTANT_NOISE, INITIAL_CONDITIONS, filename, print_start, print_interval) {
    .Call('_sRACIPEv03_multiGeneCircuit_EM_uniform_Darray_multiprint', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, D_max, D_shot_scaling, GENE_NOISE_SCALING, file_writing_interval, D_levels, D_scaling, output_precision, ANNEALING, CONSTANT_NOISE, INITIAL_CONDITIONS, filename, print_start, print_interval)
}

multiGeneCircuit_RK_adaptive_deterministic <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, RK_TOLERANCE, filename) {
    .Call('_sRACIPEv03_multiGeneCircuit_RK_adaptive_deterministic', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, RK_TOLERANCE, filename)
}

multiGeneCircuit_RK_deterministic_knockout <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, filename, knockout = NA_integer_) {
    .Call('_sRACIPEv03_multiGeneCircuit_RK_deterministic_knockout', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, filename, knockout)
}

multiGeneCircuit_RK_deterministic_multiprint <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, filename, print_start, print_interval) {
    .Call('_sRACIPEv03_multiGeneCircuit_RK_deterministic_multiprint', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, filename, print_start, print_interval)
}

multiGeneCircuit_RK_deterministic_thMod <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, filename) {
    .Call('_sRACIPEv03_multiGeneCircuit_RK_deterministic_thMod', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, filename)
}

multiGeneCircuit_RK_deterministic <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, filename) {
    .Call('_sRACIPEv03_multiGeneCircuit_RK_deterministic', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, tot_time, median_range, standard_deviation_factor, number_gene, output_precision, INITIAL_CONDITIONS, filename)
}

threshold_calculator_uniform <- function(gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, median_range, standard_deviation_factor) {
    .Call('_sRACIPEv03_threshold_calculator_uniform', PACKAGE = 'sRACIPEv03', gene_interaction, threshold_gene, g_min, g_max, k_min, k_max, possible_interactions, model_count_max, threshold_max, h, lambda_min, lambda_max, n_min, n_max, median_range, standard_deviation_factor)
}

