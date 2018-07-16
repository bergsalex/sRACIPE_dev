rm(list = ls())
library(sRACIPEv03)
working_directory <- getwd()
config_file = "inputs/sRACIPE.cfg"
configuration <- sRACIPE_load_configuration(config_file = config_file)

l1.tpos <- seq(1,1)
for(l1.tpo.counter in 1:length(l1.tpos))
{
  #l1.tpo.counter <- 1 # test
  tpo.filename <- paste("g3_",l1.tpos[l1.tpo.counter],"_1_1.txt", sep = "")
  topology.file <- file.path(getwd(),"inputs",tpo.filename)

  f.tpo <- read.table(topology.file, header = TRUE,stringsAsFactors = FALSE)
  links <- dim(f.tpo)[1]
  inhib.file <- 2^links

  for(inhib.file.counter in 1: inhib.file)
  {
    for(self.links.counter in 1:27){
    topology_file <- paste("g3_",l1.tpos[l1.tpo.counter],"_",inhib.file.counter,"_", self.links.counter,".txt", sep = "")
    print("Configuration file")
    print(topology_file)
topology_file <- paste(getwd(), "/inputs/",topology_file,sep="")
    topology <- sRACIPE_load_topology(topology_file = topology_file)

    results_directory <- ifelse(!dir.exists(file.path(working_directory, "results")), dir.create(file.path(working_directory, "results")), file.path(working_directory, "results"))

    output_file <- sRACIPE_RK_deterministic_multiprint(topology_file = topology_file, NUM_MODELS = 10, INITIAL_CONDITIONS = 1, SIM_TIME = 100.0, STEP_SIZE = 0.01, OUTPUT_PRECISION = 8,PRINT_START = 50.0, PRINT_INTERVAL = 10)

    rm(output_file)
}
  }

}









#output_file <- sRACIPE_RK_adaptive_deterministic(topology_file = topology_file)

#output_file <- sRACIPE_stochastic(topology_file = topology_file)

#output_file <- sRACIPE_RK_deterministic_knockout(topology_file = topology_file)

#output_file <- sRACIPE_RK_deterministic_knockout(topology_file = topology_file, KNOCKOUT = "GCM1",NUM_MODELS = 2000)


#plot_data_deterministic(output_file, plot_filename = topology$filename)

#plot_data_stochastic(output_file, plot_filename=topology$filename, topology_df=topology, config = configuration, bin_count=40)

#plot_interactive_heatmap(output_file, plot_filename = topology$filename)
# ###########################################################################
# #Knockout Simulations
#
# topology_file <- "/Users/koharv/Documents/Work/ScenicMay/IPSC/output/IPSC_network_w001_corr002_Target.txt"
# topology <- sRACIPE_load_topology(topology_file)
# sRACIPEv03::plot_network(topology)
# output_file <- sRACIPEv03::sRACIPE_RK_deterministic(topology_file =  topology_file, NUM_MODELS = 5000)
#
# KO = "MAZ"
# output_file_knockout <- sRACIPEv03::sRACIPE_RK_deterministic_knockout(topology_file =  topology_file, NUM_MODELS = 5000, KNOCKOUT = KO)
# sRACIPEv03::plot_data_knockout_single(output_file, output_file_knockout, KNOCKOUT = KO,plot_filename = topology$filename, topology_df = topology)
#
# #KO = "HIF1A"
# KNOCKOUT = "HMGA2"
# output_file <- "/Users/koharv/Documents/GitHub/sRACIPE_dev/results/sRACIPE_RK_IPSC_network_w001_corr002_Target_g31_output.txt"
# output_file_knockout <- paste("/Users/koharv/Documents/GitHub/sRACIPE_dev/results/sRACIPE_RK_IPSC_network_w001_corr002_Target_",KNOCKOUT,"_KO_g31_output.txt", sep = "")
# plot_data_knockout_single(output_file, output_file_knockout, KNOCKOUT = KNOCKOUT,plot_filename = topology$filename, topology_df = topology)
# #sRACIPEv03::plot_interactive_heatmap(output_file,plot_filename = topology$filename, topology_df = topology)
#
# sRACIPEv03::plot_data_deterministic(output_file,plot_filename = topology$filename, topology_df = topology)
#
# data_simulation <- read.table(output_file, header = F)
# data_simulation <- load_data(data_simulation, topology)
#
# heatmap(t(data_simulation), col=plot_color, hclustfun = function(x) hclust(x,method = 'ward.D2'), distfun=function(x) as.dist((1-cor(t(x), method = "spear"))/2))

