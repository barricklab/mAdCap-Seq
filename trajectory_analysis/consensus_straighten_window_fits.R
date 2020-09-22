#fits windows across all populations using the same fitness increase
consensus_straighten_window_fits <- function(output_path, base_file_names, fitness.bootstrap = F)
{
  #debug
  #output_path="grouped_output"
  #base_file_names= c("A1_population_window", "A2_population_window", "A3_population_window", "A6_population_window", "A7_population_window", "A9_population_window")
    
  source("common.R")
  
  ### Fitness bootstrapping
  num.bootstraps = 200
  
  ### First generation of fitting fitness
  first.generation.of.fitness.fitting = 183
  
  ################################################
  
  #Setup for graphs
  plot.width = 7
  plot.height = 5
  theme_set(theme_bw(base_size = 24))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
  line_thickness = 0.8
  log_frequency_bounds = c(-5,0)
  generation_bounds = c(163,243)
  myColors <- c("purple", "magenta", "orange", "green", "red", "blue", "brown", "cyan", "grey")
  names(myColors) <- c("hslU", "iclR", "pykF", "nadR", "spoT", "topA", "ybaL", "fabR", "other")
  colScale <- scale_colour_manual(name = "Gene",values = myColors)
  
  ### Here we infer population fitness as the best parameters that straighten out the graphs
  
  temp_merged_little_table = data.frame()
  
  for (name in base_file_names) {
    little_table = read.csv(file.path(output_path, "output", paste(name, "_significant_mutations_table.csv", sep="")))
    little_table$population = name
    temp_merged_little_table = rbind(temp_merged_little_table, little_table)
  }
  
  #set up individual fit columns
  the.generations = as.numeric(levels(as.factor(temp_merged_little_table$generation)))
  num_gen_of_data = length(the.generations)
  
  
  merged_little_table = data.frame()
  
  for (name in base_file_names) {
    little_table = read.csv(file.path(output_path, "output", paste(name, "_significant_mutations_table.csv", sep="")))
    little_table$population = name
    
    this.generations = as.numeric(levels(as.factor(little_table$generation)))
    this.num_gen_of_data = length(this.generations)
    
    for (i in 1:num_gen_of_data) {
      
      colname = paste("interval_", i, sep="")
      vector_replace = rep(0, this.num_gen_of_data)
      
      gen_final = min(the.generations[i], this.generations[length(this.generations)])
      if (i==1) {
        gen_initial = 0
      } else {
        gen_initial = the.generations[i-1]
      }
      
      k = 1
      while ((k < length(this.generations)) && (this.generations[k] < the.generations[i])) {
        k = k + 1
      }
      for (j in k:this.num_gen_of_data) {
        
        vector_replace[j] = gen_final - gen_initial
      }
      little_table[, colname] = vector_replace
    }
    
    merged_little_table = rbind(merged_little_table, little_table)
    
  }
  
  #Need to divide up mutations in different 
  merged_little_table$full_name = paste0(merged_little_table$population, "-", merged_little_table$full_name)
  
  
  write.csv(merged_little_table, file=file.path(output_path, "output", paste("consensus_significant_mutations_table.csv", sep="")))
  
  formula_string = "cbind(variant_read_count, not_variant_read_count) ~ 0 + full_name:generation + full_name"
  
 
  #################################
  
  ## Set up the real table
  little_table = merged_little_table
  
  
  
  #set up individual fit columns
  the.generations = as.numeric(levels(as.factor(little_table$generation)))
  num_gen_of_data = length(levels(as.factor(little_table$generation)))
  
  
  set_up_model_data <- function(test.population.fitness.increments, num.intervals, data.table) {
    
    #print(num.intervals)
    num.intervals.with.fitness = length(test.population.fitness.increments)
    #print(num.intervals.with.fitness)
    
    new.data.table = data.table
    
    #The effective generation is each time interval divided by the population fitness in that interval
    new.data.table[, "effective_generations"] = 0
    new.data.table[, "generation_offset"] = 0
    
    for (i in 1:num.intervals) {
      colname = paste("interval_", i, sep="")
      #print(colname)
      
      if (i>num.intervals-num.intervals.with.fitness) {
        this.fitness = 1 + sum(test.population.fitness.increments[1:(i - (num.intervals-num.intervals.with.fitness))])
      } else {
        this.fitness = 1
      }
      
      #print(i)
      #print(this.fitness)
      
      new.data.table[, "effective_generations"] = new.data.table[, "effective_generations"] + new.data.table[, colname] / this.fitness
      new.data.table[, "generation_offset"] = new.data.table[, "generation_offset"] + new.data.table[, colname] * this.fitness
    }
    
    return(new.data.table)
  }
  
  population_fitness_model <- function(test.population.fitness.increments, num.intervals, data.table) {
    
    the.data.table = set_up_model_data(test.population.fitness.increments, num.intervals, data.table)
    formula_string = "cbind(variant_read_count, not_variant_read_count) ~ 0 + full_name:effective_generations + full_name + offset(-generation)"
    fit.model = glm(as.formula(formula_string), data=the.data.table, binomial)
    return(fit.model)
  }
  
  population_fitness_model_AIC <- function(test.population.fitness.increments, num.intervals, data.table) {
    
    fit.model = population_fitness_model(test.population.fitness.increments, num.intervals, data.table)
    the.AIC = AIC(fit.model)
    return(the.AIC)
  }
  
  ## Pick mutations with resampling => note, we must give them unique names!
  bootstrap_table <- function(in.table) {
    
    mut.names = unique(in.table$full_name)
    num.muts = length(mut.names)
    this.mut.selections = sample(mut.names, length(mut.names), replace=T)
    
    new.table = data.frame()
    for (mut.name in this.mut.selections) {
      #print(mut.name)
      
      new.table = new.table %>% bind_rows(in.table %>% filter(full_name==mut.name))
    }
    
    return(new.table)
  }
  

  
  ### Bootstrap the mutations
  num.fitness.increments = sum(the.generations>first.generation.of.fitness.fitting)
  
  set_up_model_data(rep(0,num.fitness.increments), num_gen_of_data, little_table)
  
  ## Find ML values for real data
  real.data.results = optim(rep(0,num.fitness.increments), population_fitness_model_AIC, lower=rep(0,num.fitness.increments), upper=rep(0.1,num.fitness.increments), method="L-BFGS-B", num.intervals = num_gen_of_data, data.table=little_table) 
  
  real.data.fitness.increments = real.data.results$par
  real.data.fit.model = population_fitness_model(real.data.fitness.increments, num.intervals = num_gen_of_data, data.table=little_table)
  
  bootstrap.fitness.increments = matrix(rep(NA,num.fitness.increments*num.bootstraps), nrow=num.bootstraps)
  
  if (fitness.bootstrap) {
    
    ##We need to add the population to the mutation name when bootstrapping
    little_table_for_bootstrap = little_table

    for(b in 1:num.bootstraps)
    {
      cat("Bootstrap: ", b, "\n")
      this.bootstrap.table = bootstrap_table(little_table_for_bootstrap)
      
      bootstrap.results = optim(real.data.fitness.increments, population_fitness_model_AIC, lower=rep(0,num.fitness.increments), upper=rep(0.1,num.fitness.increments), method="L-BFGS-B", num.intervals = num_gen_of_data, data.table=this.bootstrap.table) 
      
      print(bootstrap.results$par)
      bootstrap.fitness.increments[b,] = bootstrap.results$par
    }
  } else {
    for(b in 1:num.bootstraps) {
      bootstrap.fitness.increments[b,] = real.data.fitness.increments
    }
  }
  
  #Now add up to actual fitness values
  bootstrap.fitness = bootstrap.fitness.increments
  
  for (c in 1:ncol(bootstrap.fitness)) {
    for (r in 1:nrow(bootstrap.fitness)) {
      bootstrap.fitness[r,c] = sum(bootstrap.fitness.increments[r,1:c])
    }
  }
  
  bootstrap.fitness.CI95 = data.frame()
  for (c in 1:ncol(bootstrap.fitness)) {
    q = quantile(bootstrap.fitness[,c], probs=c(0.025, 0.975))
    
    bootstrap.fitness.CI95 = rbind(bootstrap.fitness.CI95,
                                   data.frame(generation.end = the.generations[num_gen_of_data-num.fitness.increments+c], 
                                              population.relative.fitness.L95 = q[1], 
                                              population.relative.fitness.U95 =q[2])
    )
    
  }
  
  #Add maximum likelihood values
  bootstrap.fitness.CI95$population.relative.fitness = NA
  for (r in 1:length(real.data.fitness.increments)) {
    bootstrap.fitness.CI95$population.relative.fitness[r] = sum(real.data.fitness.increments[1:r])
  }
  
  # Convert to relative fitness
  bootstrap.fitness.CI95$population.relative.fitness = 1+bootstrap.fitness.CI95$population.relative.fitness/log(2)
  bootstrap.fitness.CI95$population.relative.fitness.L95 = 1+bootstrap.fitness.CI95$population.relative.fitness.L95/log(2)
  bootstrap.fitness.CI95$population.relative.fitness.U95 = 1+bootstrap.fitness.CI95$population.relative.fitness.U95/log(2)
  rownames(bootstrap.fitness.CI95) <- NULL
  
  #Initial rows
  
  initial.rows = data.frame(
    generation.end=the.generations[1:(num_gen_of_data-num.fitness.increments)],
    population.relative.fitness=rep(1,num_gen_of_data-num.fitness.increments),
    population.relative.fitness.L95=rep(NA,num_gen_of_data-num.fitness.increments),
    population.relative.fitness.U95=rep(NA,num_gen_of_data-num.fitness.increments)
  )
  
  bootstrap.fitness.CI95 = initial.rows %>% bind_rows(bootstrap.fitness.CI95)
  bootstrap.fitness.CI95$generation.start = 0
  bootstrap.fitness.CI95$generation.start[2:length(bootstrap.fitness.CI95$generation.start)] = bootstrap.fitness.CI95$generation.end[1:(length(bootstrap.fitness.CI95$generation.start)-1)]
  bootstrap.fitness.CI95 = bootstrap.fitness.CI95 %>% relocate(generation.start, .before = generation.end)
  
  bootstrap.fitness.CI95$was.fit = (bootstrap.fitness.CI95$generation.end > the.generations[num_gen_of_data-num.fitness.increments]) 
  
  write.csv(bootstrap.fitness.CI95, file=file.path(output_path, "output", paste("consensus_population_fitness.csv", sep="")), row.names=F)
  
  
  ## Step plot of the fitness over time
  
  graphing.fitness.step = bootstrap.fitness.CI95
  graphing.fitness.step = graphing.fitness.step %>% rbind(graphing.fitness.step[nrow(graphing.fitness.step),])
  graphing.fitness.step$generation.start[nrow(graphing.fitness.step)] = graphing.fitness.step$generation.end[nrow(graphing.fitness.step)]
  graphing.fitness.step$population.relative.fitness.L95[nrow(graphing.fitness.step)] = NA
  graphing.fitness.step$population.relative.fitness.U95[nrow(graphing.fitness.step)] = NA
  graphing.fitness.step$generation.mean = (graphing.fitness.step$generation.start + graphing.fitness.step$generation.end)/2
  graphing.fitness.step$generation = graphing.fitness.step$generation.start
  
  p = ggplot(graphing.fitness.step) + 
    geom_step(aes(x=generation, y=population.relative.fitness), color="black", direction = 'hv') + 
    geom_errorbar(aes(x=generation.mean, ymin=population.relative.fitness.L95, ymax=population.relative.fitness.U95), width=2) + 
    coord_cartesian(xlim=c(the.generations[1],the.generations[num_gen_of_data]), ylim=c(1,1.12))
  ggsave(filename=file.path(output_path, "output", paste("consensus_population_fitness_step.pdf", sep="")), width=plot.width, height=plot.height)
  
  #####################
  
  ## compare to per population fitness model with same number of parameters
  
#  formula_string = "cbind(variant_read_count, not_variant_read_count) ~ 0 + full_name:generation + full_name"
#  for (i in best_model_index:num_gen_of_data) {
#    formula_string = paste(formula_string, " + population:fitness_", as.character(i), sep="")
#  }
#  per_population_model = glm(as.formula(formula_string), data=merged_little_table, binomial)
  
#  anova(best_model, per_population_model, test="Chisq")
  
  ##Per-population models fit much better
  
  best_model = real.data.fit.model
  best_model.confint = confint(best_model)
  set_up_model_table = set_up_model_data(real.data.fitness.increments, num_gen_of_data, little_table)
  
  
  population.accounted.for.fitness.table = data.frame()
  
  for (base_file_name in base_file_names) {
    
    this.population.table = subset(merged_little_table, population==base_file_name)
    
    straightened.filtered.output.table = read.csv(file=file.path(output_path, "output", paste(base_file_name, "_significant_mutations.csv", sep="")))
    straightened.filtered.output.table$original.selection.coefficient = straightened.filtered.output.table$selection.coefficient
    straightened.filtered.output.table$original.selection.coefficient.CI95L = straightened.filtered.output.table$selection.coefficient.CI95L
    straightened.filtered.output.table$original.selection.coefficient.CI95U = straightened.filtered.output.table$selection.coefficient.CI95U
    
    straightened.filtered.output.table$original.slope = straightened.filtered.output.table$slope
    straightened.filtered.output.table$original.slope.stderr = straightened.filtered.output.table$slope.stderr
    
    for (i in 1:nrow(straightened.filtered.output.table)) {
      
      this_full_name = as.character(straightened.filtered.output.table$full_name[i])
      this_full_name_plus_population = paste0(base_file_name, "-" , this_full_name)
      
      slope_name = paste("full_name", this_full_name_plus_population, ":effective_generations", sep="")
      intercept_name = paste("full_name", this_full_name_plus_population, sep="")
      
      this.slope = coef(best_model)[slope_name]
      this.slope.stderr =  summary(best_model)$coefficients[, 2][slope_name]
      
      straightened.filtered.output.table$fitness.effect[i] = (this.slope-1)/log(2)
      straightened.filtered.output.table$fitness.effect.CI95L[i] = (best_model.confint[slope_name, 1]-1)/log(2)
      straightened.filtered.output.table$fitness.effect.CI95U[i] = (best_model.confint[slope_name, 2]-1)/log(2)
      
      ##graph the straightened fits
      # this section produces graphs for all of the ones that were significant that include
      # binomial confidence intervals that account for the counts
      test_series = subset(set_up_model_table, full_name==this_full_name_plus_population)
      
      test_series$fit = predict(best_model, test_series, type="response")
      test_series$log_fit = log10(predict(best_model, test_series, type="response"))
      
      p = ggplot(test_series, aes(x=generation, y=log_frequency, color=gene, group=evidence))
      
      for(i in 1:nrow(test_series)) {
        bt = binom.test(test_series$variant_read_count[i], test_series$total_read_count[i])
        test_series$L95.CI[i] = log10(bt$conf.int[1] * test_series$freq.corr[i])
        test_series$U95.CI[i] = log10(bt$conf.int[2] * test_series$freq.corr[i])
      }
      
      if (all(test_series$variant_read_count != test_series$total_read_count)) {    
        p = ggplot(test_series, aes(x=generation, y=log_frequency, color=gene, group=evidence))
        p + geom_line(size=line_thickness) + geom_line(aes(x=generation, y=log_fit), color="black") + geom_errorbar(aes(ymin=L95.CI, ymax=U95.CI), width=0.2) + coord_cartesian(xlim=generation_bounds, ylim=log_frequency_bounds)       
        ggsave(filename=file.path(output_path, "graphs", paste(base_file_name, "/", this_full_name, "_consensus.straightened", ".pdf", sep="")), width=plot.width, height=plot.height)
      }
      
    }
    
    write.csv(straightened.filtered.output.table, file=file.path(output_path, "output", paste(base_file_name, "_consensus_straightened_significant_mutations.csv", sep="")), row.names=F)
    
    #Now calculate the theoretical population fitness from the observed variants and what proportion it explains
    name_to_fitness = straightened.filtered.output.table %>% select(full_name, fitness.effect)
    merged_little_table_plus_fitness =  merged_little_table %>% filter(population==base_file_name) %>% left_join(name_to_fitness, by="full_name")
    
    for (i in 1:(nrow(graphing.fitness.step)-1)) {
      g = graphing.fitness.step$generation.end[i]
      
      this_gen = merged_little_table_plus_fitness %>% filter(generation==g)

      if (nrow(this_gen)==0) {
        next
      }
      
      population.accounted.for.fitness.table = rbind(population.accounted.for.fitness.table,
                           data.frame(
                             population = base_file_name, 
                             generation=g,
                             relative.fitness = sum(this_gen$frequency*this_gen$fitness.effect)
                           )
      )
    }
  }
  
  write.csv(population.accounted.for.fitness.table, file.path(output_path, "output", paste("consensus_tracked_mutation_total_fitness.csv", sep="")))
}

