rare_graphs_complete <- function (variant_read_counts_file_name, total_read_counts_file_name, output_path, base_file_name, autocorrelation_threshold)
{
  
  #autocorrelation_threshold = 0.55
  
  source("common.R")
  
  #Setup for graphs
  plot.width = 7
  plot.height = 5
  theme_set(theme_bw(base_size = 24))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
  line_thickness = 0.8
  log_frequency_bounds = c(-4,0)
  expanded_log_frequency_bounds = c(-5,0)
  generation_bounds = c(30,537)
  myColors <- c("green", "orange", "blue", "red", "cyan", "brown", "purple", "magenta", "grey")
  names(myColors) <- c("nadR", "pykF", "topA", "spoT", "fabR", "ybaL", "iclR", "hslU", "other")
  colScale <- scale_colour_manual(name = "Gene",values = myColors)
  
  dir.create(file.path(output_path), showWarnings = FALSE)
  
  dir.create(file.path(output_path, "output"), showWarnings = FALSE)
  dir.create(file.path(output_path, "graphs"), showWarnings = FALSE)
  dir.create(file.path(output_path, "graphs", base_file_name), showWarnings = FALSE)
  
  ### also allow us to only look at certain generational timepoints (remove the ones that were collected in a different Illumina run or failed)
  remove_generations = c()
  if (base_file_name=="A1_population_complete") {
    remove_generations = c(169,183,196,209,223,236)
  } else if (base_file_name=="A2_population_complete") {
    remove_generations = c(169,183,196,209,223,236,270.12,510.24)
  } else  if (base_file_name=="A3_population_complete") {
    remove_generations = c(169,183,196,209,223,236)
  } else  if (base_file_name=="A7_population_complete") {
    remove_generations = c(169,183,196,209,223,236)
  }
  
  res = load_count_files(variant_read_counts_file_name, total_read_counts_file_name)
  mutation_info <- res$mutation_info 
  big_table <- res$big_table
  
  #Remove certain generations where data is missing or from different illumina runs
  big_table = subset(big_table, !(generation %in% remove_generations))
  
  min_log_value = -6
  big_table$log_frequency[big_table$log_frequency==-Inf] <- min_log_value
  
  p = ggplot(big_table, aes(x=generation, y=log_frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + colScale + coord_cartesian(xlim=generation_bounds, ylim=log_frequency_bounds)
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_all_mutations_log_frequency.pdf", sep="")), width=plot.width, height=plot.height)
  
  p = ggplot(big_table, aes(x=generation, y=frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + colScale + coord_cartesian(xlim=generation_bounds, ylim=c(0, 1))
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_all_mutations_frequency.pdf", sep="")), width=plot.width, height=plot.height)
  
  write.csv(big_table, file=file.path(output_path, "output", paste(base_file_name, "_full_table.csv", sep="")), row.names=F)
  
  
  ### Look for autocorrelation here
  
  full.output.table = data.frame()
  required_nonzero_observations = 2
  filtered_table = data.frame();
  filtered_table_2 = data.frame();
  passing_full_name_indices = c()
  for (i in 1:nrow(mutation_info)) {
    test_series = subset(big_table, full_name==mutation_info$full_name[i])
    
    the.autocorrelation.coefficient = NA
    
    #require at least three nonzero values
    nonzero_series = subset(test_series, variant_read_count>0)
    the.non_zero_observations = nrow(nonzero_series)
    
    if (the.non_zero_observations >= required_nonzero_observations)  {
      
      the.autocorrelation.coefficient = autocorrelation_coefficient(test_series$frequency)
      
      if ((!is.na(the.autocorrelation.coefficient)) && (the.autocorrelation.coefficient >= autocorrelation_threshold) ) {
        filtered_table = rbind(filtered_table, test_series)
        passing_full_name_indices = c(passing_full_name_indices, i)
      }
    }
    
    new_row = data.frame(full_name=mutation_info$full_name[i], gene=mutation_info$gene[i], type=mutation_info$type[i], category=mutation_info$category[i], position=mutation_info$position[i], non.zero = the.non_zero_observations, autocorrelation.coefficient=the.autocorrelation.coefficient)
    full.output.table = rbind (full.output.table, new_row)
  }
  write.csv(full.output.table, file=file.path(output_path, "output", paste(base_file_name, "_all_mutations.csv", sep="")), row.names=F)
  
  
  total_table = data.frame()
  for (i in levels(as.factor(filtered_table$generation))) {
    total_row = data.frame(generation=as.numeric(i), frequency = sum((subset(filtered_table,generation==i))$frequency) )
    total_table = rbind(total_table, total_row)
  }
  total_table$log_frequency = log10(total_table$frequency)
  total_table$evidence = 1;
  
  write.csv(total_table, file=file.path(output_path, "output", paste(base_file_name, "_mutation_totals.csv", sep="")), row.names=F)
  
  nrow(filtered_table)
  p = ggplot(filtered_table, aes(x=generation, y=log_frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + colScale + coord_cartesian(xlim=generation_bounds, ylim=log_frequency_bounds)
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_log_frequency.pdf", sep="")), width=plot.width, height=plot.height)
  
  p = ggplot(filtered_table, aes(x=generation, y=log_frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + geom_line(data=total_table, aes(x=generation, y=log_frequency), size=line_thickness, color="black") + colScale + coord_cartesian(xlim=generation_bounds, ylim=log_frequency_bounds)
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_log_frequency_including_total.pdf", sep="")), width=plot.width, height=plot.height)
  
  p = ggplot(filtered_table, aes(x=generation, y=frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + colScale + coord_cartesian(xlim=generation_bounds, ylim=c(0, 1))
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_frequency.pdf", sep="")), width=plot.width, height=plot.height)
  
  p = ggplot(filtered_table, aes(x=generation, y=frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + geom_line(data=total_table, aes(x=generation, y=frequency), size=line_thickness, color="black") + colScale + coord_cartesian(xlim=generation_bounds, ylim=c(0, 1))
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_frequency_including_total.pdf", sep="")), width=plot.width, height=plot.height)
  
  write.csv(filtered_table, file=file.path(output_path, "output", paste(base_file_name, "_filtered_mutations.csv", sep="")), row.names=F)
  
  output.table = data.frame()
  for(n in 1:length(passing_full_name_indices)) {
    mutation_info_index = passing_full_name_indices[n]
    this_full_name = as.character(mutation_info$full_name[passing_full_name_indices[n]])
    # this section produces graphs for all of the ones that were significant that include
    # binomial confidence intervals that account for the counts
    test_series = subset(filtered_table, full_name==this_full_name)
    
    the.autocorrelation.coefficient = autocorrelation_coefficient(test_series$frequency)
    
    for(i in 1:nrow(test_series)) {
      bt = binom.test(test_series$variant_read_count[i], test_series$total_read_count[i])
      test_series$L95.CI[i] = log10(bt$conf.int[1] * test_series$freq.corr[i]) 
      test_series$U95.CI[i] = log10(bt$conf.int[2] * test_series$freq.corr[i]) 
    }
    
    if (all(test_series$variant_read_count != test_series$total_read_count)) {    
      p = ggplot(test_series, aes(x=generation, y=log_frequency, color=gene, group=evidence))
      p + geom_line(size=line_thickness) + geom_errorbar(aes(ymin=L95.CI, ymax=U95.CI), width=0.2) + coord_cartesian(xlim=generation_bounds, ylim=expanded_log_frequency_bounds)
      ggsave(filename=file.path(output_path, "graphs", paste(base_file_name, "/", this_full_name, ".pdf", sep="")), width=plot.width, height=plot.height) 
    }
    
    new_row = data.frame(full_name=mutation_info$full_name[mutation_info_index], gene=mutation_info$gene[mutation_info_index], type=mutation_info$type[mutation_info_index], category=mutation_info$category[mutation_info_index], position=mutation_info$position[mutation_info_index], autocorrelation.coefficient=the.autocorrelation.coefficient, original.index=mutation_info_index)
    output.table = rbind (output.table, new_row)
  }
  
  write.csv(output.table, file=file.path(output_path, "output", paste(base_file_name, "_significant_mutations.csv", sep="")), row.names=F)
  
  
}
