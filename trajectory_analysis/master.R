## Setup for run
source("common.R")
source("rare_graphs_complete.R")
source("rare_graphs_window.R")
source("grouped_mutations_merge.R")
source("grouped_mutations_find.R")
source("consensus_straighten_window_fits.R")
source("detection_limits.R")
source("combination_graphs.R")


complete_base_file_names = c("A1_population_complete", "A2_population_complete", "A3_population_complete", "A7_population_complete")

window_base_file_names = c("A1_population_window", "A2_population_window", "A3_population_window", "A6_population_window", "A7_population_window", "A9_population_window")

input_path = "initial_input"
output_path = "initial_output"

grouped_input_path = "grouped_input"
grouped_output_path = "grouped_output"

this_autocorrelation_threshold = 0.55

for (base_file_name in complete_base_file_names) {
  cat(base_file_name, "\n")
  rare_graphs_complete(file.path(input_path, paste(base_file_name, ".variant.tsv", sep="")), file.path(input_path, paste(base_file_name, ".total.tsv", sep="")), output_path, base_file_name, this_autocorrelation_threshold)
}

for (base_file_name in window_base_file_names) {
  cat(base_file_name, "\n")
  rare_graphs_window(file.path(input_path, paste(base_file_name, ".variant.tsv", sep="")), file.path(input_path, paste(base_file_name, ".total.tsv", sep="")), output_path, base_file_name)
}


for (base_file_name in c(window_base_file_names, complete_base_file_names)) {
  cat(base_file_name, "\n")
  grouped_mutations_find(file.path(input_path, paste(base_file_name, ".variant.tsv", sep="")), file.path(input_path, paste(base_file_name, ".total.tsv", sep="")), output_path, base_file_name)
}


#####
# Manual curation step
#
# Copy the files 'initial_output/output/<base_name>_grouped_mutations.csv' and edit to group mutations to meet desired criteria
# Store the curated files as 'grouped_mutations/<base_name>_grouped_mutations.csv'
# The scripts in group_mutations_helper_workspace are used for creating valid mutation entries from this evidence
# Graphs of the grouped evidence are generated under grouped_output/graphs/grouping.
#####


# Combine the variant/total counts for mutations that are designated as merged
grouped_mutations_path = "grouped_mutations"
for (base_file_name in c(window_base_file_names, complete_base_file_names)) {
  
  small_base_name = base_file_name
  small_base_name = gsub("_window", "", small_base_name)
  small_base_name = gsub("_complete", "", small_base_name)
  
  cat(small_base_name, "\n")
  
  grouped_mutations_merge(input_path, paste(base_file_name, ".variant.tsv", sep=""), file.path(grouped_mutations_path, paste(small_base_name, "_grouped_mutations.csv", sep="")), grouped_input_path)
  grouped_mutations_merge(input_path, paste(base_file_name, ".total.tsv", sep=""), file.path(grouped_mutations_path, paste(small_base_name, "_grouped_mutations.csv", sep="")), grouped_input_path)
}


this_autocorrelation_threshold = 0.55
for (base_file_name in complete_base_file_names) {
    cat(base_file_name, "\n")
    rare_graphs_complete(file.path(grouped_input_path, paste(base_file_name, ".variant.tsv", sep="")), file.path(grouped_input_path, paste(base_file_name, ".total.tsv", sep="")), paste(grouped_output_path, sep=""), base_file_name, this_autocorrelation_threshold)
}

for (base_file_name in window_base_file_names) {
  cat(base_file_name, "\n")
  rare_graphs_window(file.path(grouped_input_path, paste(base_file_name, ".variant.tsv", sep="")), file.path(grouped_input_path, paste(base_file_name, ".total.tsv", sep="")), grouped_output_path, base_file_name, fitness.bootstrap=T)
}

consensus_straighten_window_fits(grouped_output_path, window_base_file_names)

combination_graphs()

detection_limits()
