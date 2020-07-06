source("common.R")
library(dplyr)
options(dplyr.summarise.inform=F) 

detection_limits_samples <- function(samples) {

  cutoff = 0.05
  complete.table = data.frame()
  complete.gen.above.freq.table = data.frame()
  full.sig.table = data.frame()
  
  for (base.name in samples) {
    full.table = read.csv(file.path("grouped_output", "output", paste0(base.name, "_full_table.csv")))
    sig.muts = read.csv(file.path("grouped_output", "output", paste0(base.name, "_significant_mutations.csv")))
    sig.table = full.table %>% filter(full_name %in% unique(sig.muts$full_name))

    this.population = base.name
    this.population = sub("_complete", "", this.population)
    this.population = sub("_window", "", this.population)
    
    sig.table$population = this.population
    full.sig.table = full.sig.table %>% bind_rows(sig.table)
    
    max.freq.table = sig.table %>% group_by(full_name) %>% summarize(max.frequency = max(frequency))
    max.freq.table$population = this.population
    
    complete.table = complete.table %>% bind_rows(max.freq.table)
    gen.above.freq.table = sig.table %>% filter(frequency >= 0.05) %>% group_by(full_name) %>% summarize(min.generation = min(generation), max.generation = max(generation))
    gen.above.freq.table$gen.above = gen.above.freq.table$max.generation - gen.above.freq.table$min.generation
    gen.above.freq.table$population = this.population
    
    complete.gen.above.freq.table = complete.gen.above.freq.table %>% bind_rows(gen.above.freq.table)
  }

  print(complete.gen.above.freq.table)
  
  complete.table$log10.max.frequency = log10(complete.table$max.frequency )
  ggplot(complete.table, aes(log10.max.frequency)) + stat_ecdf(geom = "step")
  
  complete.frequency.table = data.frame()
  thresholds = c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1)
  for (threshold in thresholds) {
    this.frequency.table = complete.table %>% filter(max.frequency >= threshold) %>% group_by(full_name, population) %>% summarize(max.frequency = max(max.frequency)) %>% ungroup() %>% summarize(n = n())
    this.frequency.table$threshold = threshold
    complete.frequency.table = complete.frequency.table %>% bind_rows(this.frequency.table)
  }
  
  complete.frequency.table$never.reached = max(complete.frequency.table$n) - complete.frequency.table$n
  
  print(complete.frequency.table)
  
  full.sig.table %>% filter(full_name == "REL606_nadr.RA.2280.A.G.0.gene_name.nadR.snp_type.nonsynonymous") %>% filter(frequency > 0.0) %>% group_by(full_name, population) %>% summarize(first.detected.generation = min(generation))
  
  ##  For those that surpass X% how much earlier did we detect at all?
  thresholds = c(0, 0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1)
  for (threshold in thresholds) {
    this.first.detected.table = full.sig.table %>% filter(frequency > 0.0) %>% group_by(full_name, population) %>% summarize(first.detected.generation = min(generation))
    this.first.above.threshold.table = full.sig.table %>% filter(frequency >= threshold) %>% group_by(full_name, population) %>% summarize(first.above.threshold.generation = min(generation))
    
    threshold.table = this.first.above.threshold.table %>% left_join(this.first.detected.table, by=c("full_name", "population"))
    
    mean.generations.earlier = mean(threshold.table$first.above.threshold.generation - threshold.table$first.detected.generation)
    cat("Mean number of generations detected before threshold: ", threshold, " = ", mean.generations.earlier, "\n")
  }

}

detection_limits <- function()
{
  
  complete_base_file_names = c("A1_population_complete", "A2_population_complete", "A3_population_complete", "A7_population_complete")
  
  window_base_file_names = c("A1_population_window", "A2_population_window", "A3_population_window", "A6_population_window", "A7_population_window", "A9_population_window")
  
  
  cat("Complete time course samples\n")
  
  detection_limits_samples(complete_base_file_names)
  
  cat("Window time course samples\n")
  
  detection_limits_samples(window_base_file_names)
  
  
  cat("All samples\n")
  
  detection_limits_samples(c(complete_base_file_names, window_base_file_names))

  
}
