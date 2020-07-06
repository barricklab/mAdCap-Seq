#fits windows across all populations using the same fitness increase
consensus_straighten_window_fits <- function(output_path, base_file_names)
{
  source("common.R")
  
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
  
  
  output_path = grouped_output_path
  base_file_names = window_base_file_names
  
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
    
    for (i in 2:num_gen_of_data) {
      
      colname = paste("fitness_", i, sep="")
      vector_replace = rep(0, this.num_gen_of_data)
      
      gen_final = min(the.generations[i], this.generations[length(this.generations)])
      gen_initial = the.generations[i-1]
      
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
  
  write.csv(merged_little_table, file=file.path(output_path, "output", paste("consensus_significant_mutations_table.csv", sep="")))
  
  formula_string = "cbind(variant_read_count, not_variant_read_count) ~ 0 + full_name:generation + full_name"
  
  best_model_index = 0
  best.AIC = 99999999999
  best_total_fitnesses = c()
  fit_slope_intercept = c()
  for (i in 1:num_gen_of_data) {
    
    fit_slope_intercept[[i]] = glm(as.formula(formula_string), data=merged_little_table, binomial)
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
  best_model.confint = confint(best_model, c("fitness_6", "fitness_7", "fitness_8", "fitness_9", "fitness_10"))
  
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
    ggsave(filename=file.path(output_path, "output", paste("consensus_population_fitness_step.pdf", sep="")), width=plot.width, height=plot.height)
    
    #bar plot version
    p = ggplot(fitnesses)
    p + geom_bar(aes(x=generation, y=fitness, width = w), position = "identity", stat = "identity", fill="gray", color="black") + geom_errorbar(aes(x=generation, ymin=L, ymax=U), width=2) + coord_cartesian(xlim=c(the.generations[1],the.generations[num_gen_of_data]), ylim=c(-0.1,0.15))
    ggsave(filename=file.path(output_path, "output", paste("consensus_population_fitness.pdf", sep="")), width=plot.width, height=plot.height)
    
  }
  
  ## compare to per population fitness model with same number of parameters
  
  formula_string = "cbind(variant_read_count, not_variant_read_count) ~ 0 + full_name:generation + full_name"
  for (i in best_model_index:num_gen_of_data) {
    formula_string = paste(formula_string, " + population:fitness_", as.character(i), sep="")
  }
  per_population_model = glm(as.formula(formula_string), data=merged_little_table, binomial)
  
  anova(best_model, per_population_model, test="Chisq")
  
  ##Per-population models fit much better
  
  
  for (base_file_name in base_file_names) {
    
    this.population.table = subset(merged_little_table, population==base_file_name)
    
    best_model.confint = confint(best_model, paste0("full_name", as.factor(this.population.table$full_name), ":generation") )
    
    straightened.filtered.output.table = read.csv(file=file.path(output_path, "output", paste(base_file_name, "_significant_mutations.csv", sep="")))
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
      test_series = subset(this.population.table, full_name==this_full_name)
      
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
        ggsave(filename=file.path(output_path, "graphs", paste(base_file_name, "/", this_full_name, ".consensus.straightened", ".pdf", sep="")), width=plot.width, height=plot.height)
      }
      
    }
    
    write.csv(straightened.filtered.output.table, file=file.path(output_path, "output", paste(base_file_name, "_consensus_straightened_significant_mutations.csv", sep="")), row.names=F)
    
  }
}