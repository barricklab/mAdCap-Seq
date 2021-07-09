source("common.R")
library(dplyr)
options(dplyr.summarise.inform=F) 

fraction_fitness_increase <- function() {
  
  window_base_file_names = c("A1_population_window", "A2_population_window", "A3_population_window", "A6_population_window", "A7_population_window", "A9_population_window")
  
  
  for (generation.of.interest in c(196, 209, 223)) {
    
    increases = c();
    
    for (base.name in window_base_file_names) {
      
      full.table = read.csv(file.path("grouped_output", "output", paste0(base.name, "_full_table.csv")))
      sig.muts = read.csv(file.path("grouped_output", "output", paste0(base.name, "_consensus_straightened_significant_mutations.csv")))
      
      pop.fitness = read.csv(file.path("grouped_output", "output", paste0(base.name, "_population_fitness.csv")))
      
      # Take the population fitness at this time point as the average of the windows before and after
      pop.fitness.increase = ( (pop.fitness %>% filter(generation.start == generation.of.interest))$population.relative.fitness[1] +
                                 (pop.fitness %>% filter(generation.end == generation.of.interest))$population.relative.fitness[1] ) / 2 - 1
      
      #Filter to just the ones with estimated fitness at the generation of interest
      sig.table = full.table %>% 
        filter(generation == generation.of.interest) %>% 
        filter(full_name %in% unique(sig.muts$full_name))
      
      #Join in the fitness values
      sig.table = sig.table %>% left_join(sig.muts, by="full_name")
      
      #Calculate the total fitness we accounted for
      mutation.fitness.increase = sum(sig.table$fitness.effect * sig.table$frequency)
      
      fitness.increase.percent = round(mutation.fitness.increase / pop.fitness.increase, 4) * 100
      cat(base.name, "\n")
      cat("  ", "Generation                 : ", generation.of.interest, "\n")
      cat("  ", "Num mutations tracked      : ", nrow(sig.table), "\n")
      cat("  ", "Population fitness increase: ", pop.fitness.increase, "\n")
      cat("  ", "Mutation fitness increase  : ", mutation.fitness.increase, "\n")
      cat("  ", "Percent accounted for      : ", fitness.increase.percent, "%\n")
      
      increases = c(increases,fitness.increase.percent)
    }
    
    cat("AVERAGE INCREASE  : ", mean(increases), "%\n")
    cat("MIN INCREASE      : ", min(increases), "%\n")
    cat("MAX INCREASE      : ", max(increases), "%\n")
  }
 
  
  base_file_names = c("A1_population", "A2_population", "A3_population", "A7_population")
  generations.of.interest = c(296.80)
  pop.fitness.increase = 0.09
  
  increases = c();
  
  
  for (base.name in base_file_names) {
    
    window_base_file_name = paste0(base.name, "_window")
    complete_base_file_name = paste0(base.name, "_complete")
    
    sig.muts = read.csv(file.path("grouped_output", "output", paste0(window_base_file_name, "_consensus_straightened_significant_mutations.csv")))
    full.table = read.csv(file.path("grouped_output", "output", paste0(complete_base_file_name, "_full_table.csv")))

    complete.sig.muts = read.csv(file.path("grouped_output", "output", paste0(complete_base_file_name, "_significant_mutations.csv")))
    
    
    full.table = full.table %>% filter(generation %in% generations.of.interest) 
    
    sig.table = full.table %>%filter(full_name %in% unique(sig.muts$full_name))
    
    sig.table = sig.table %>% left_join(sig.muts, by="full_name")
    
    mutation.fitness.increase = sum(sig.table$fitness.effect * sig.table$frequency)
    
    fitness.increase.percent = round(mutation.fitness.increase / pop.fitness.increase, 4) * 100
    
    all.sig.table = full.table %>%filter(full_name %in% unique(complete.sig.muts$full_name))
    
    all.mutations.percent = round(sum(all.sig.table$frequency), 4) * 100 
    
    
    cat(base.name, "\n")
    cat("  ", "Generation                 : ", generations.of.interest, "\n")
    cat("  ", "Num mutations tracked      : ", nrow(sig.table), "\n")
    cat("  ", "Population fitness increase: ", pop.fitness.increase, "\n")
    cat("  ", "Mutation fitness increase  : ", mutation.fitness.increase, "\n")
    cat("  ", "Percent accounted for      : ", fitness.increase.percent, "%\n")
    cat("  ", "All mutations percent      : ", all.mutations.percent, "%\n")
    
    increases = c(increases,fitness.increase.percent)
    
  }
  
  cat("AVERAGE INCREASE  : ", mean(increases), "%\n")
  cat("MIN INCREASE      : ", min(increases), "%\n")
  cat("MAX INCREASE      : ", max(increases), "%\n")
   
}
