rare_graphs_window <- function (variant_read_counts_file_name, total_read_counts_file_name, output_path, base_file_name, fitness.bootstrap=F)
{
  source("common.R")
  
  ## For testing, set these
  #variant_read_counts_file_name = "initial_input/A6_population_window.variant.tsv"
  #total_read_counts_file_name = "initial_input/A6_population_window.total.tsv"
  #output_path = "initial_output"
  #base_file_name = "A6_population_window"
  ### End testing setup
  
  # DM 100 is 2E8 cells/ml, the culture size is 10 ml, and there are 1:100 dilutions
  effective.population.size = 1/((1/2E7 + 1/2E9)/2)
  
  ########### Filter settings ##################
  
  ### Taken into account before fitting slopes
  minimum_nonzero_observations_required = 0
  
  final_frequency_required = 0.0001
  final_frequency_required_at_generation_and_later = 223
  
  minimum_first_nonzero_observation_generation = 196
  
  minimum_y_intercept = -20
  
  ### Also required to build graphs
  slope_p_value_required = 0.005
  
  ### Also required to be retained as "significant"
  maximum_AIC_required = 200

  ### Fitness bootstrapping
  num.bootstraps = 1000
  
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
  generation_bounds = c(163,243.44)
  myColors <- c("purple", "magenta", "orange", "green", "red", "blue", "brown", "cyan", "grey")
  names(myColors) <- c("hslU", "iclR", "pykF", "nadR", "spoT", "topA", "ybaL", "fabR", "other")
  colScale <- scale_colour_manual(name = "Gene",values = myColors)
  
  
  dir.create(file.path(output_path), showWarnings = FALSE)
  
  dir.create(file.path(output_path, "graphs"), showWarnings = FALSE)
  dir.create(file.path(output_path, "output"), showWarnings = FALSE)
  dir.create(file.path(output_path, "graphs", base_file_name), showWarnings = FALSE)
  
  res = load_count_files(variant_read_counts_file_name, total_read_counts_file_name)
  mutation_info <- res$mutation_info 
  big_table <- res$big_table
  num.time.points = nrow(big_table)/nrow(mutation_info)
  
  p = ggplot(big_table, aes(x=generation, y=log_frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + colScale + coord_cartesian(xlim=generation_bounds, ylim=log_frequency_bounds)
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_all_mutations_log_frequency.pdf", sep="")), width=plot.width, height=plot.height)
  
  write.csv(big_table, file=file.path(output_path, "output", paste(base_file_name, "_full_table.csv", sep="")), row.names=F)
  
  ### first rule out entries based on not having enough observations or not a sufficient frequency
  ### so that these are not included in the p-value adjustment
  
  new_mutation_info = data.frame() 
  for (i in 1:nrow(mutation_info)) {
    test_series = subset(big_table, full_name==mutation_info$full_name[i])
    
    test_series = test_series
    
    #number of nonzero points
    nonzero_series = subset(test_series, variant_read_count>0)
    
    #first observation
    passed_generation_test =  min(nonzero_series$generation) <= minimum_first_nonzero_observation_generation
    the.non_zero_observations = nrow(nonzero_series)
    
    #observations past a given generation must have a certain frequency
    not.meeting.test = test_series %>% filter(generation >= final_frequency_required_at_generation_and_later) %>% filter(is.na(frequency)  || (frequency < final_frequency_required) )
    passed_final_frequency_test = nrow(not.meeting.test) == 0
        
    #cat(failed_final_frequency_test, "\n")
    
    if ((the.non_zero_observations >= minimum_nonzero_observations_required) && passed_final_frequency_test && passed_generation_test) {
      new_mutation_info = rbind(new_mutation_info, mutation_info[i,])
    }
  }
  
  cat("Requiring enough non-zero observation and final frequency cutoff:\n")
  cat("  Candidate mutations reduced from ", nrow(mutation_info), " to ", nrow(new_mutation_info), "\n")
  mutation_info = new_mutation_info
  
  full.output.table = data.frame()
  filtered.output.table = data.frame()
  filtered_table = data.frame();
  filtered_table_2 = data.frame();
  passing_full_name_indices = c()
  for (i in 1:nrow(mutation_info)) {
    
    test_series = subset(big_table, full_name==mutation_info$full_name[i])
    the.p_value = NA
    the.adjusted.p_value = NA
    the.intercept = NA
    the.intercept.stderr = NA
    the.slope = NA
    the.slope.stderr = NA
    the.deviance = NA
    the.AIC = NA
    
    fit_slope_intercept = glm(cbind(variant_read_count, not_variant_read_count) ~ generation + offset(-generation), binomial, data=test_series)
    
    fit_intercept = glm(cbind(variant_read_count, not_variant_read_count) ~ 1, binomial, data=test_series)
    
    significant_slope_anova = anova(fit_intercept,fit_slope_intercept,test="Chi")
    
    the.intercept = coef(fit_slope_intercept)[1]
    the.intercept.stderr = (summary(fit_slope_intercept))$coefficients[1,2]
    the.slope = coef(fit_slope_intercept)[2]
    the.slope.stderr = (summary(fit_slope_intercept))$coefficients[2,2]
    
    the.fitness.effect = (the.slope-1)/log(2)
    #the.fitness.effect.stderr = (summary(fit_slope_intercept))$coefficients[2,2]/log(2)
    
    the.AIC = fit_slope_intercept$aic
    the.deviance = fit_slope_intercept$deviance
    
    the.p_value = significant_slope_anova$"Pr(>Chi)"[2]

    ### Extra info
    the.generation.first.appearance = (-log(effective.population.size) - the.intercept)/the.slope
    nonzero_series = subset(test_series, variant_read_count>0)
    the.non_zero_observations = nrow(nonzero_series)
    
    # Make row to output
    
    new_row = data.frame(full_name=mutation_info$full_name[i], gene=mutation_info$gene[i], type=mutation_info$type[i], category=mutation_info$category[i], position=mutation_info$position[i], non.zero = the.non_zero_observations, intercept=the.intercept, intercept.stderr=the.intercept.stderr, slope=the.slope, slope.stderr=the.slope.stderr, slope.p_value=the.p_value, fitness.effect=the.fitness.effect, deviance=the.deviance, AIC=the.AIC, generation.first.appearance=the.generation.first.appearance)
    full.output.table = rbind (full.output.table, new_row)

    if ((the.AIC <= maximum_AIC_required) && (the.fitness.effect>0) && (the.p_value<=slope_p_value_required) && (the.intercept >= minimum_y_intercept)) {
      filtered_table_2 = rbind(filtered_table_2, test_series)
      filtered.output.table = rbind (filtered.output.table, new_row)          
    }
    
  }
  
  cat("Requiring fitness effect > 0 and AIC cutoff <= ", maximum_AIC_required, " and p-value for nonzero slope of <= ", slope_p_value_required, " and minimum y-intercept > ", minimum_y_intercept, ":\n")
  cat("  Candidate mutations reduced from ", nrow(new_mutation_info), " to ", nrow(filtered.output.table), "\n")

  #### Could add filtering on adjusted P.value here
  #filtered.output.table$slope.adjusted.p_value = p.adjust(filtered.output.table$slope.p_value, method="BH")
  
  ## write full lists with info
  write.csv(full.output.table, file=file.path(output_path, "output", paste(base_file_name, "_all_mutations.csv", sep="")), row.names=F)
  
  write.csv(filtered.output.table, file=file.path(output_path, "output", paste(base_file_name, "_significant_mutations.csv", sep="")), row.names=F)
  
  ## graph overall distributions of what we are filtering on
  ggplot(full.output.table, aes(x=AIC)) + geom_histogram() + coord_cartesian()  
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_AIC_distribution.pdf", sep=""))) 
  
  ggplot(full.output.table, aes(x=intercept)) + geom_histogram() + coord_cartesian()  
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_intercept_distribution.pdf", sep=""))) 
  
  ggplot(full.output.table, aes(x=slope)) + geom_histogram() + coord_cartesian()  
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_slope_distribution.pdf", sep=""))) 
  
  
  ## graph the fits
  nrow(filtered_table_2)
  p = ggplot(filtered_table_2, aes(x=generation, y=log_frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + colScale + coord_cartesian(xlim=generation_bounds, ylim=log_frequency_bounds)
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_log_frequency.pdf", sep="")), width=plot.width, height=plot.height) 
  
  total_table = data.frame()
  for (i in levels(as.factor(filtered_table_2$generation))) {
    total_row = data.frame(generation=as.numeric(i), frequency = sum((subset(filtered_table_2,generation==i))$frequency) )
    total_table = rbind(total_table, total_row)
  }
  total_table$log_frequency = log10(total_table$frequency)
  total_table$evidence = 1;

  p = ggplot(rbind(filtered_table_2), aes(x=generation, y=log_frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + geom_line(data=total_table, aes(x=generation, y=log_frequency), size=line_thickness, color="black") + colScale + coord_cartesian(xlim=generation_bounds, ylim=log_frequency_bounds)
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_log_frequency_including_total.pdf", sep="")), width=plot.width, height=plot.height) 
  
  nrow(filtered_table_2)
  p = ggplot(filtered_table_2, aes(x=generation, y=frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + colScale + coord_cartesian(xlim=generation_bounds, ylim=c(0,1))
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_frequency.pdf", sep="")), width=plot.width, height=plot.height)
  
  nrow(filtered_table_2)
  p = ggplot(filtered_table_2, aes(x=generation, y=frequency, color=gene, group=evidence))
  p + geom_line(size=line_thickness) + geom_line(data=total_table, aes(x=generation, y=frequency), size=line_thickness, color="black") + colScale + coord_cartesian(xlim=generation_bounds, ylim=c(0,1))
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_frequency_including_total.pdf", sep="")), width=plot.width, height=plot.height) 
  
  output.table = data.frame()
  for(this_full_name in unique(filtered_table_2$full_name)) {
    
    # this section produces graphs for all of the ones that were significant that include
    # binomial confidence intervals that account for the counts
    test_series = subset(filtered_table_2, full_name==this_full_name)
    
    fit_slope_intercept = glm(cbind(test_series$variant_read_count, test_series$not_variant_read_count) ~ test_series$generation, binomial)
    fit_intercept = glm(cbind(test_series$variant_read_count, test_series$not_variant_read_count) ~ 1, binomial)
    significant_slope_anova = anova(fit_intercept,fit_slope_intercept,test="Chi")
    the.p_value = significant_slope_anova$"Pr(>Chi)"[2]
    the.adjusted.p_value = the.p_value * nrow(mutation_info)
    the.AIC = fit_slope_intercept$aic
    the.deviance = fit_slope_intercept$deviance
    the.intercept = coef(fit_slope_intercept)[1]
    the.intercept.stderr = (summary(fit_slope_intercept))$coefficients[1,2]
    the.slope = coef(fit_slope_intercept)[2]
    the.slope.stderr = (summary(fit_slope_intercept))$coefficients[2,2]
    
    for(i in 1:nrow(test_series)) {
      bt = binom.test(test_series$variant_read_count[i], test_series$total_read_count[i])
      test_series$L95.CI[i] = log10(bt$conf.int[1] * test_series$freq.corr[i])
      test_series$U95.CI[i] = log10(bt$conf.int[2] * test_series$freq.corr[i] )
    }
    
    test_series$log_fit = log10(predict(fit_slope_intercept, test_series, type="response"))
    test_series$fit = predict(fit_slope_intercept, test_series, type="response")
    
    if (all(test_series$variant_read_count != test_series$total_read_count)) {    
      p = ggplot(test_series, aes(x=generation, y=log_frequency, color=gene, group=evidence))
      p + geom_line(size=line_thickness) + geom_line(aes(x=generation, y=log_fit), color="black") + geom_errorbar(aes(ymin=L95.CI, ymax=U95.CI), width=0.2) + coord_cartesian(xlim=generation_bounds, ylim=log_frequency_bounds)       
      ggsave(filename=file.path(output_path, "graphs", paste(base_file_name, "/", this_full_name, ".pdf", sep="")), width=plot.width, height=plot.height)
    }
  }
  
  
  write.csv(filtered_table_2, file=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_table.csv", sep="")))

  
  ### Here we infer population fitness as the best parameters that straighten out the graphs
  
  ######################
  
  #Fitness must increase over time, so all are constrained to be positive
  
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
  
  ## Set up the real table
  little_table = filtered_table_2
  
  #set up individual fit columns
  the.generations = as.numeric(levels(as.factor(little_table$generation)))
  num_gen_of_data = length(levels(as.factor(little_table$generation)))
  
  little_table[, "interval_1"] = the.generations[1] 
  
  for (i in 2:num_gen_of_data) {
    colname = paste("interval_", i, sep="")
    vector_replace = rep(0, num_gen_of_data)
    
    gen_final = the.generations[i]
    gen_initial = the.generations[i-1]
    for (j in i:num_gen_of_data) {
      vector_replace[j] = gen_final - gen_initial
    }
    little_table[, colname] = vector_replace
  }
    
  
  ### Bootstrap the mutations
  num.fitness.increments = sum(the.generations>first.generation.of.fitness.fitting)
  
  ## Find ML values for real data
  real.data.results = optim(rep(0,num.fitness.increments), population_fitness_model_AIC, lower=rep(0,num.fitness.increments), upper=rep(0.1,num.fitness.increments), method="L-BFGS-B", num.intervals = num_gen_of_data, data.table=little_table) 
  
  real.data.fitness.increments = real.data.results$par
  real.data.fit.model = population_fitness_model(real.data.fitness.increments, num.intervals = num_gen_of_data, data.table=little_table)
  
  bootstrap.fitness.increments = matrix(rep(NA,num.fitness.increments*num.bootstraps), nrow=num.bootstraps)
  
  if (fitness.bootstrap) {
    
    for(b in 1:num.bootstraps)
    {
      cat("Bootstrap: ", b, "\n")
      this.bootstrap.table = bootstrap_table(little_table)
      
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
  
  write.csv(bootstrap.fitness.CI95, file=file.path(output_path, "output", paste(base_file_name, "_population_fitness.csv", sep="")), row.names=F)
  
  
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
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_population_fitness_step.pdf", sep="")), width=plot.width, height=plot.height)
  

  ## Graph the straightened fit output
  
  best_model = real.data.fit.model
  best_model.confint = confint(best_model)
  
  straightened.filtered.output.table = filtered.output.table
  straightened.filtered.output.table$original.fitness.effect = straightened.filtered.output.table$fitness.effect

  straightened.filtered.output.table$original.slope = straightened.filtered.output.table$slope
  straightened.filtered.output.table$original.slope.stderr = straightened.filtered.output.table$slope.stderr
  
  set_up_model_table = set_up_model_data(real.data.fitness.increments, num_gen_of_data, little_table)
  
  for (i in 1:nrow(straightened.filtered.output.table)) {
    
    this_full_name = as.character(straightened.filtered.output.table$full_name[i])
    
    slope_name = paste("full_name", this_full_name, ":effective_generations", sep="")
    intercept_name = paste("full_name", this_full_name, sep="")
    
    this.slope = coef(best_model)[slope_name]
    this.slope.stderr =  summary(best_model)$coefficients[, 2][slope_name]
    
    straightened.filtered.output.table$fitness.effect[i] = (this.slope-1)/log(2)
    straightened.filtered.output.table$fitness.effect.CI95L[i] = (best_model.confint[slope_name, 1]-1)/log(2)
    straightened.filtered.output.table$fitness.effect.CI95U[i] = (best_model.confint[slope_name, 2]-1)/log(2)
    
    ##graph the straightened fits
    # this section produces graphs for all of the ones that were significant that include
    # binomial confidence intervals that account for the counts
    test_series = subset(set_up_model_table, full_name==this_full_name)
    
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
      ggsave(filename=file.path(output_path, "graphs", paste(base_file_name, "/", this_full_name, ".straightened", ".pdf", sep="")), width=plot.width, height=plot.height)
    }
  
  }
  
  write.csv(straightened.filtered.output.table, file=file.path(output_path, "output", paste(base_file_name, "_straightened_significant_mutations.csv", sep="")), row.names=F)
  
}


