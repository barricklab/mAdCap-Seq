rare_graphs_window <- function (variant_read_counts_file_name, total_read_counts_file_name, output_path, base_file_name)
{
  source("common.R")
  
  ## For testing, set these
  #variant_read_counts_file_name = "grouped_input/A1_population_window.variant.tsv"
  #total_read_counts_file_name = "grouped_input/A1_population_window.total.tsv"
  #output_path = "grouped_output"
  #base_file_name = "A1_population_window"
  ### End testing setup
  
  # DM 100 should be 2E8, but pop size is smaller due to bottlenecks
  population.size = 2E6
  
  ########### Filter settings ##################
  
  ### Taken into account before fitting slopes and not used when correcting p-values
  minimum_nonzero_observations_required = 2
  final_frequency_required = 0.0001
  final_frequency_points_required = 3
  
  ### Also required to build graphs
  slope_p_value_required = 0.05
  
  ### Also required to be retained as "significant"
  maximum_AIC_required = 200
  minimum_intercept_required = -15
  
  ################################################
  
  #Setup for graphs
  plot.width = 7
  plot.height = 5
  theme_set(theme_bw(base_size = 24))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
  line_thickness = 0.8
  log_frequency_bounds = c(-5,0)
  generation_bounds = c(133,213)
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
  
  ### first rule out entries base on not having enough observations
  ### so that these are not included in the p-value adjustment
  
  new_mutation_info = data.frame() 
  for (i in 1:nrow(mutation_info)) {
    test_series = subset(big_table, full_name==mutation_info$full_name[i])
    
    nonzero_series = subset(test_series, variant_read_count>0)
    the.non_zero_observations = nrow(nonzero_series)
    
    #require at least a certain number of the last timepoints to be above some threshold frequency
    failed_final_frequency_test = 0
    
    test_generations = test_series$generation
    test_generations = sort(test_generations, decreasing=T)
    
    for (u in 1:min(length(test_generations), final_frequency_points_required)) {
      test_row = subset(test_series, generation == test_generations[u])
      #cat(as.character(test_row$full_name[1]), " ", test_generations[u], " ", test_row$frequency, "\n")
      if ( (is.na(test_row$frequency[1])) || (test_row$frequency[1] < final_frequency_required) ) {
        failed_final_frequency_test = 1
      }
    }
    #cat(failed_final_frequency_test, "\n")
    
    if ((the.non_zero_observations >= minimum_nonzero_observations_required) && (failed_final_frequency_test == 0)) {
      new_mutation_info = rbind(new_mutation_info, mutation_info[i,])
    }
  }
  
  cat("rows changed from ", nrow(mutation_info), " to ", nrow(new_mutation_info), "\n")
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
    
    
    ##### Unused tests
    #require first data point be below 1%
    first_data_point = test_series$log_frequency[1]
    
    #special code for ruling out first run artifact
    g1_log_freq = subset(test_series, generation==133.4)$log_frequency
    g2_log_freq = subset(test_series, generation==160.08)$log_frequency
    g3_log_freq = subset(test_series, generation==186.76)$log_frequency
    g4_log_freq = subset(test_series, generation==213.44)$log_frequency
    
    min_first_run_requirement_test = min(g1_log_freq, g2_log_freq, g3_log_freq, g4_log_freq)
    #if (is.na(first_data_point) || (first_data_point < -2)) {
    ####################
    
    
    fit_slope_intercept = glm(cbind(test_series$variant_read_count, test_series$not_variant_read_count) ~ test_series$generation, binomial)
    
    fit_intercept = glm(cbind(test_series$variant_read_count, test_series$not_variant_read_count) ~ 1, binomial)
    
    significant_slope_anova = anova(fit_intercept,fit_slope_intercept,test="Chi")
    
    the.intercept = coef(fit_slope_intercept)[1]
    the.intercept.stderr = (summary(fit_slope_intercept))$coefficients[1,2]
    the.slope = coef(fit_slope_intercept)[2]
    the.slope.stderr = (summary(fit_slope_intercept))$coefficients[2,2]
    the.AIC = fit_slope_intercept$aic
    the.deviance = fit_slope_intercept$deviance
    
    # Bonferroni corrected
    the.p_value = significant_slope_anova$"Pr(>Chi)"[2]
    the.adjusted.p_value = the.p_value * nrow(mutation_info)
    
    ### Extra info
    the.generation.first.appearance = (-log(population.size) - the.intercept)/the.slope
    nonzero_series = subset(test_series, variant_read_count>0)
    the.non_zero_observations = nrow(nonzero_series)
    
    # Make row to output
    
    new_row = data.frame(full_name=mutation_info$full_name[i], gene=mutation_info$gene[i], type=mutation_info$type[i], category=mutation_info$category[i], position=mutation_info$position[i], non.zero = the.non_zero_observations, intercept=the.intercept, intercept.stderr=the.intercept.stderr, slope=the.slope, slope.stderr=the.slope.stderr,  selection.coefficient=exp(the.slope)-1, slope.p_value=the.p_value, slope.adjusted.p_value=the.adjusted.p_value, deviance=the.deviance, AIC=the.AIC, generation.first.appearance=the.generation.first.appearance)
    full.output.table = rbind (full.output.table, new_row)
    
    #the.adjusted.p_value
    if (the.adjusted.p_value  <= slope_p_value_required) {
      filtered_table = rbind(filtered_table, test_series)
      passing_full_name_indices = c(passing_full_name_indices, i)
      
      if ((the.intercept >= minimum_intercept_required) && (the.AIC <= maximum_AIC_required)) {
        filtered_table_2 = rbind(filtered_table_2, test_series)
        filtered.output.table = rbind (filtered.output.table, new_row)          
      }
    }
    
  }
  
  ## write full lists with info
  write.csv(full.output.table, file=file.path(output_path, "output", paste(base_file_name, "_all_mutations.csv", sep="")), row.names=F)
  
  write.csv(filtered.output.table, file=file.path(output_path, "output", paste(base_file_name, "_significant_mutations.csv", sep="")), row.names=F)
  
  ## graph overall distributions of what we are filtering on
  ggplot(subset(full.output.table, slope.adjusted.p_value <= slope_p_value_required), aes(x=AIC)) + geom_histogram() + coord_cartesian()  
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_AIC_distribution.pdf", sep=""))) 
  
  ggplot(subset(full.output.table, slope.adjusted.p_value <= slope_p_value_required), aes(x=intercept)) + geom_histogram() + coord_cartesian()  
  ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_intercept_distribution.pdf", sep=""))) 
  
  ggplot(subset(full.output.table, slope.adjusted.p_value <= slope_p_value_required), aes(x=slope)) + geom_histogram() + coord_cartesian()  
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
  for(n in 1:length(passing_full_name_indices)) {
    mutation_info_index = passing_full_name_indices[n]
    this_full_name = as.character(mutation_info$full_name[passing_full_name_indices[n]])
    # this section produces graphs for all of the ones that were significant that include
    # binomial confidence intervals that account for the counts
    test_series = subset(filtered_table, full_name==this_full_name)
    
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
  
### Here we infer population fitness as the best parameters that straighten out the graphs
  
  write.csv(filtered_table_2, file=file.path(output_path, "output", paste(base_file_name, "_significant_mutations_table.csv", sep="")))

  little_table = filtered_table_2
  
  #set up individual fit columns
  the.generations = as.numeric(levels(as.factor(little_table$generation)))
  num_gen_of_data = length(levels(as.factor(little_table$generation)))
  for (i in 2:num_gen_of_data) {
    colname = paste("fitness_", i, sep="")
    vector_replace = rep(0, num_gen_of_data)
    
    gen_final = the.generations[i]
    gen_initial = the.generations[i-1]
    for (j in i:num_gen_of_data) {
      vector_replace[j] = gen_final - gen_initial
    }
    little_table[, colname] = vector_replace
  }
  
  formula_string = "cbind(variant_read_count, not_variant_read_count) ~ 0 + full_name:generation + full_name"
  
  best_model_index = 0
  best.AIC = 99999999999
  best_total_fitnesses = c()
  fit_slope_intercept = c()
  for (i in 1:num_gen_of_data) {
    
    fit_slope_intercept[[i]] = glm(as.formula(formula_string), data=little_table, binomial)
    #a = anova(fit_slope_intercept[[i]], fit_slope_intercept[[i+1]], test = "Chisq")
    cat("Model with ", i-1, " fitness parameters ::: \n");
    
    the.AIC = AIC(fit_slope_intercept[[i]] )
    K = attr(logLik(fit_slope_intercept[[i]] ), "df")
    n = nobs(fit_slope_intercept[[i]] )
    the.AICc = the.AIC + 2 * K * (K + 1) / (n - K - 1)
    
    cat ("  AIC = ", the.AIC, " AICc = ", the.AICc, "\n")
    
    first_fit_gen = the.generations[num_gen_of_data-i+1]
    
    all_fitnesses = rep(0, num_gen_of_data-1)
    if (i != 1) {

      
      
      for (j in (i-1):1) {
        parameter_name = paste("fitness_", as.character(num_gen_of_data-j+1), sep="")
        gen_final = the.generations[num_gen_of_data-j + 1]
        gen_initial = the.generations[num_gen_of_data-j]
        all_fitnesses[num_gen_of_data-j] =  -(coef(fit_slope_intercept[[i]])[parameter_name])
        
        cat("  ", parameter_name, " = ", coef(fit_slope_intercept[[i]])[parameter_name], " gen (", gen_initial, "-", gen_final, ") total fitness in interval (per gen) = ", -(coef(fit_slope_intercept[[i]])[parameter_name]), "\n")
      }
    }
    
    
    #if ((the.AIC < best.AIC) && (sum(all_fitnesses) > 0)) {
    #if ((the.AIC < best.AIC)) {
    if (first_fit_gen == 166)
    {
      cat("   ====> Best because first_fit_gen == 166", "\n")
      
      best_model_index = i
      best.AIC = the.AIC;
      #best_total_fitnesses = total_fitnesses;
    }
    
    formula_string = paste(formula_string, " + fitness_", as.character(num_gen_of_data-i+1), sep="")
  }
  
  ## now we need to extract the parameters from the best model, save, and graph
  
  cat("Best model has this many fitness parameters: ", (best_model_index-1), "\n")

  best_model = fit_slope_intercept[[best_model_index]]
  best_model.confint = confint(best_model)

  ## here's a bar graph of the fitness over time with error bars

  
  if (best_model_index != 1) {
    
    fitnesses = data.frame()
    step_fitnesses = data.frame()
    step_errors = data.frame()
    
    step_fitnesses = rbind(step_fitnesses, data.frame(generation=the.generations[1], fitness=0))
    
    this.fitness = 0
    for (i in (num_gen_of_data-best_model_index+2):num_gen_of_data ) {
    
      colname = paste("fitness_", i, sep="")

      this.fitness = exp(-coef(best_model)[colname])-1
      this.fitness.confint.U = exp(-best_model.confint[colname,1])-1
      this.fitness.confint.L = exp(-best_model.confint[colname,2])-1
      
      final_gen = the.generations[i]
      initial_gen = the.generations[i-1]
      midpoint_gen = (final_gen + initial_gen) / 2
      
      fitnesses = rbind(fitnesses, data.frame(generation=midpoint_gen, w = final_gen-initial_gen, fitness=this.fitness, U = this.fitness.confint.U, L = this.fitness.confint.L))
      
      step_fitnesses = rbind(step_fitnesses, data.frame(generation=initial_gen, fitness=this.fitness))
      
      step_errors = rbind(step_errors, data.frame(generation=midpoint_gen, L = this.fitness.confint.L, U = this.fitness.confint.U))
    }
    step_fitnesses = rbind(step_fitnesses, data.frame(generation=the.generations[num_gen_of_data], fitness=this.fitness))

    #step plot version
    
    p = ggplot(step_errors)
    p + geom_step(data=step_fitnesses, aes(x=generation, y=fitness), color="black", direction = 'hv') + geom_errorbar(aes(x=generation, ymin=L, ymax=U), width=2) + coord_cartesian(xlim=c(the.generations[1],the.generations[num_gen_of_data]), ylim=c(-0.1,0.15))
    ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_population_fitness_step.pdf", sep="")), width=plot.width, height=plot.height)
    
    #bar plot version
    p = ggplot(fitnesses)
    p + geom_bar(aes(x=generation, y=fitness, width = w), position = "identity", stat = "identity", fill="gray", color="black") + geom_errorbar(aes(x=generation, ymin=L, ymax=U), width=2) + coord_cartesian(xlim=c(the.generations[1],the.generations[num_gen_of_data]), ylim=c(-0.1,0.15))
    ggsave(filename=file.path(output_path, "output", paste(base_file_name, "_population_fitness.pdf", sep="")), width=plot.width, height=plot.height)
    
  }
  
  straightened.filtered.output.table = filtered.output.table
  straightened.filtered.output.table$original.selection.coefficient = straightened.filtered.output.table$selection.coefficient
  straightened.filtered.output.table$selection.coefficient.CI95L = straightened.filtered.output.table$selection.coefficient.stderr
  straightened.filtered.output.table$selection.coefficient.CI95U = straightened.filtered.output.table$selection.coefficient.stderr
  
  straightened.filtered.output.table$original.slope = straightened.filtered.output.table$slope
  straightened.filtered.output.table$original.slope.stderr = straightened.filtered.output.table$slope.stderr
  
  for (i in 1:nrow(straightened.filtered.output.table)) {
    
    this_full_name = as.character(straightened.filtered.output.table$full_name[i])
    
    slope_name = paste("full_name", this_full_name, ":generation", sep="")
    intercept_name = paste("full_name", this_full_name, sep="")
    
    this.slope = coef(best_model)[slope_name]
    this.slope.stderr =  summary(best_model)$coefficients[, 2][slope_name]
    
    straightened.filtered.output.table$selection.coefficient[i] = exp(this.slope)-1
    straightened.filtered.output.table$selection.coefficient.CI95L[i] = exp(best_model.confint[slope_name, 1])-1
    straightened.filtered.output.table$selection.coefficient.CI95U[i] = exp(best_model.confint[slope_name, 2])-1
    
    straightened.filtered.output.table$slope[i] = this.slope
    straightened.filtered.output.table$slope.stderr[i] = this.slope.stderr 
    
    ##graph the straightened fits
    # this section produces graphs for all of the ones that were significant that include
    # binomial confidence intervals that account for the counts
    test_series = subset(little_table, full_name==this_full_name)
    
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


