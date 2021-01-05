combination_graphs <- function () 
{
  
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(scales) 
  
  output_path = "grouped_output"
  summary_output_path = file.path(output_path, "summary")
  dir.create(summary_output_path, showWarnings = FALSE)
  
  ###############################################################################
  ## Global setup for graphing
  ############################################################################### 
  
  myTypeColors <- c("orange", "black", "gray", "blue", "red", "green", "purple", "magenta")
  names(myTypeColors) <- c("IS_element", "indel", "big_indel", "nonsynonymous", "nonsense", "synonymous","intergenic", "junction")
  typeColScale <- scale_colour_manual(name = "type",values = myTypeColors)
  typeColFillScale = scale_fill_manual(name = "gene",values = myTypeColors)
  
  #size in nt of the region input as template for read mapping!
  gene_sizes = c( pykF=4443, fabR=3781, ybaL=4533, hslU=4738, nadR=4032, iclR=3710, spoT=5231, topA=5656)
  gene_strands = c( pykF=+1, fabR=+1, ybaL=-1, hslU=-1, nadR=+1, iclR=-1, spoT=+1, topA=5656)
  gene_starts = c( pykF=1631, fabR=1677, ybaL=1400, hslU=1400, nadR=1400, iclR=1400, spoT=1429, topA=1659)
  gene_ends = c( pykF=3043, fabR=2381, ybaL=3076, hslU=2731, nadR=2632, iclR=2224, spoT=3831, topA=4256)

  theme_set(theme_bw(base_size = 24))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
  myColors <- c("green", "orange", "blue", "red", "cyan", "brown", "purple", "magenta", "grey")
  names(myColors) <- c("nadR", "pykF", "topA", "spoT", "fabR", "ybaL", "iclR", "hslU", "other")
  colScale <- scale_colour_manual(name = "Gene",values = myColors)
  colFillScale = scale_fill_manual(name = "Gene",values = myColors)
  
  ###############################################################################
  ## Load COMPLETE timecourse mutations
  ###############################################################################
  
  complete.file.names=c("A1_population_complete_significant_mutations.csv", "A2_population_complete_significant_mutations.csv", "A3_population_complete_significant_mutations.csv", "A7_population_complete_significant_mutations.csv")
  complete.names=c("A1", "A2", "A3", "A7")
  
  complete.table = data.frame()
  for (i in 1:length(complete.file.names)) {
    
    this.table = read.csv(file.path(output_path, "output", complete.file.names[i]), header=T)
    this.table$population = complete.names[i]
    complete.table = rbind(complete.table, this.table)
  }
  complete.table$full_name_population = paste(complete.table$full_name,complete.table$population , sep=".")
  
  complete.table$fitness.effect = NA
  complete.table$fitness.effect.CI95L = NA
  complete.table$fitness.effect.CI95U = NA
  
  write.csv(complete.table, file.path(output_path, "output", paste("complete_mutations.csv", sep="")), row.names=F)
  cat("Number of mutations identified in the complete timecourses:", nrow(complete.table))
  
  ###############################################################################
  ## Load WINDOW timecourse mutations
  ###############################################################################
  
  window.names=c("A1", "A2", "A3", "A6", "A7", "A9")
  
  window.file.names=c("A1_population_window_consensus_straightened_significant_mutations.csv", "A2_population_window_consensus_straightened_significant_mutations.csv", "A3_population_window_consensus_straightened_significant_mutations.csv", "A6_population_window_consensus_straightened_significant_mutations.csv", "A7_population_window_consensus_straightened_significant_mutations.csv", "A9_population_window_consensus_straightened_significant_mutations.csv")
  
  #Alternative analysis that does not use consensus straightened (only straightened on a per populaton basis)
  #window.file.names=c("A1_population_window_straightened_significant_mutations.csv", "A2_population_window_straightened_significant_mutations.csv", "A3_population_window_straightened_significant_mutations.csv", "A6_population_window_straightened_significant_mutations.csv", "A7_population_window_straightened_significant_mutations.csv", "A9_population_window_straightened_significant_mutations.csv")
  
  window.table = data.frame()
  for (i in 1:length(window.file.names)) {
    
    this.table = read.csv(file.path(output_path, "output", window.file.names[i]), header=T)
    this.table$population = window.names[i]
    window.table = rbind(window.table, this.table)
  }
  window.table$full_name_population = paste(window.table$full_name,window.table$population , sep=".")
  window.table$gene = factor(as.character(window.table$gene), levels=names(myColors))
  window.table$type = factor(as.character(window.table$type), levels=names(myTypeColors))
  
  write.csv(window.table, file.path(output_path, "output", paste("window_mutations.csv", sep="")), row.names=F)
  write.csv(window.table, file.path(summary_output_path, paste("window_mutations.csv", sep="")), row.names=F)
  cat("Number of mutations identified in the window timecourses:", nrow(window.table))
  
  
  ###############################################################################
  ## Create a combined.all table. For analyzing *independent* mutations 
  ## (counting multiple times if same one occurred in multiple populations!)
  ###############################################################################
  
  union.table = table(c(as.character(window.table$full_name_population), as.character(complete.table$full_name_population)))
  
  # Start table with ones only found in the window set
  window.only.table = window.table
  window.only.table$found.in.both = as.numeric(union.table[window.only.table$full_name_population])
  window.only.table = subset (window.only.table, found.in.both==1)
  window.only.table.population.totals = window.only.table %>% group_by(population) %>% summarize(n=n())
  cat("Number of mutations identified ONLY in the window timecourses:", nrow(window.only.table))
  print(window.only.table.population.totals)
  
  complete.only.table = complete.table
  complete.only.table$found.in.both = as.numeric(union.table[complete.only.table$full_name_population])
  complete.only.table = subset (complete.only.table, found.in.both==1)
  complete.only.table.population.totals = complete.only.table %>% group_by(population) %>% summarize(n=n())
  cat("Number of mutations identified ONLY in the complete timecourses:", nrow(complete.only.table))
  print(complete.only.table.population.totals)
  
  window.table = window.table %>% ungroup()
  complete.only.table = complete.only.table %>% ungroup()
  
  combined.table = data.frame(
    gene = as.factor(c(as.character(window.table$gene), as.character(complete.only.table$gene))), 
    fitness.effect = c(window.table$fitness.effect, complete.only.table$fitness.effect), 
    fitness.effect.CI95L = c(window.table$fitness.effect.CI95L, complete.only.table$fitness.effect.CI95L), 
    fitness.effect.CI95U = c(window.table$fitness.effect.CI95U, complete.only.table$fitness.effect.CI95U), 
    type = as.factor(c(as.character(window.table$type), as.character(complete.only.table$type))), 
    population = c(window.table$population, complete.only.table$population), 
    position = c(window.table$position, complete.only.table$position), 
    full_name = as.factor(c(as.character(window.table$full_name), as.character(complete.only.table$full_name))), 
    full_name_population = as.factor(c(as.character(window.table$full_name_population), as.character(complete.only.table$full_name_population))) )
  
  combined.table$found.in.both  = as.factor(union.table[combined.table$full_name_population])
  combined.table.population.totals = combined.table %>% group_by(population) %>% summarize(n=n())
  
  both.table = combined.table %>% filter(found.in.both==2)
  both.table.population.totals = both.table %>% group_by(population) %>% summarize(n=n())
  cat("Number of mutations identified in BOTH window and complete timecourses:", nrow(both.table))
  print(both.table.population.totals)
  
  combined.table.population.totals = combined.table %>% group_by(population) %>% summarize(n=n())
  cat("Number of mutations identified TOTAL in all data:", nrow(combined.table))
  print(combined.table.population.totals)
  
  combined.table.gene.totals = combined.table %>% group_by(gene) %>% summarize(n=n())
  print(combined.table.gene.totals)
  
  write.csv(combined.table, file=file.path(summary_output_path, "combined_mutations.csv"))
  
  
  
  ###############################################################################
  ## Start graphing WINDOW timecourse mutations...
  ###############################################################################
  
  ### Setting this line makes all of the following graphs for the window timecourses!
  filtered.table = window.table
  
  this.max.y = 25
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=gene)) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + colFillScale + scale_y_continuous(expand = c(0,0), limits = c(0,this.max.y))
  ggsave(file.path(summary_output_path,"window_fitness_effects_by_gene.pdf"), width=8, height=5)
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=gene)) + geom_histogram(binwidth=0.0025, alpha=1, position="stack", boundary = 0.5) + 
    colFillScale + 
    scale_y_continuous(expand = c(0,0), limits = c(0,24)) + 
    facet_grid(rows=vars(gene))
  ggsave(file.path(summary_output_path,"window_fitness_effects_by_gene.pdf"), width=8, height=5)
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=type)) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + scale_y_continuous(expand = c(0,0), limits = c(0,this.max.y)) + typeColFillScale
  ggsave(file.path(summary_output_path,"window_fitness_effects_by_type.pdf"), width=8, height=5) 
  
  facet_table = subset(filtered.table, type %in% c("IS_element", "indel", "nonsynonymous"))
  facet_table = subset(facet_table, gene %in% c("nadR", "pykF", "topA"))
  ggplot(facet_table, aes(x=fitness.effect, fill=population)) +facet_grid(type~gene) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + scale_y_continuous(expand = c(0,0)) 
  ggsave(file.path(summary_output_path,"window_selection_coefficients_by_type_by_gene.pdf"), width=8, height=5) 
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=population)) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + scale_y_continuous(expand = c(0,0))
  ggsave(file.path(summary_output_path,"window_selection_coefficients_by_population.pdf"), width=8, height=5)
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=type)) + geom_histogram(binwidth=0.002, alpha=1, position="stack", boundary=0) + scale_y_continuous(expand = c(0,0)) + facet_wrap(~population) + typeColFillScale
  ggsave(file.path(summary_output_path,"window_selection_coefficients_by_population_by_type_colors.pdf"), width=8, height=5)
  
  ggplot(filtered.table, aes(population, fill=gene)) + colFillScale + geom_bar(colour="black") + scale_y_continuous(expand = c(0,0))
  ggsave(file.path(summary_output_path,"window_selection_coefficients_fit_by_gene_by_population.pdf"), width=8, height=5)
  
  
  ####### stats interlude
  
  ## Calc average and sd of all selection coefficients
  mean(filtered.table$fitness.effect)
  sd(filtered.table$fitness.effect)
  
  
  #######
  # How many times did we see the same WINDOW mutation in different populations?
  
  mutations.in.populations = filtered.table %>% group_by(full_name) %>% summarize(num.populations=n(), populations=paste0(population, collapse = "-")) 
  
  mutations.in.multiple.populations = mutations.in.populations %>% filter(num.populations > 1)
  
  write.csv(mutations.in.multiple.populations, file.path(summary_output_path,"window_mutations_in_multiple_populations.csv"))
  
  cat("Number of window mutations that were found in more than one population: ", nrow(mutations.in.multiple.populations))
  cat("Number of window mutations (counting multi-populations multiple times): ", nrow(filtered.table))
  cat("Number of distinct window mutations (counting multi-populations once): ", nrow(filtered.table) - sum(mutations.in.multiple.populations$num.populations-1 ))
  
  
  ## Add data to the table and rewrite....
  filtered.table = filtered.table %>% left_join(mutations.in.populations, by="full_name")
  
  write.csv(filtered.table, file.path(output_path, "output", paste("window_mutations.csv", sep="")), row.names=F)
  
  
  #############
  # Is there a difference in the distribution of fitness effects across populations?
  
  kruskal.test(fitness.effect ~ population, data=filtered.table)
  
  # Yes, but it could be for stochastic reasons on the high side and
  # sampling depth on the low side. We find many more mutation in A7 
  # and those go to lower selection coefficients for example
  
  ggplot(filtered.table, aes(x=population, y=fitness.effect)) + 
    #colFillScale + 
    geom_violin(trim=FALSE) + 
    geom_jitter(width=0.2)
  
  ggsave(file.path(summary_output_path,"window_violin_fitness_effects_by_population.pdf"), width=8, height=5)
  
  ###############################################################################
  ## Stats and graphs that include only WINDOW timecourse mutations
  ###############################################################################
  
  three.genes = subset(filtered.table, ((gene == "nadR") | (gene == "pykF") | (gene=="topA")))
  ggplot( three.genes, aes(x=gene, fill=type)) + typeColFillScale + geom_bar(colour="black") + scale_y_continuous(expand = c(0,0))
  ggsave(file.path(summary_output_path,"window_mutation_types_by_gene.pdf"), width=7.5, height=5)
  
  ## Print the number of window mutations in each gene
  print(filtered.table %>% group_by(gene) %>% summarize(n=n()))
  
  ## What percent are in the top three genes
  cat("Fraction of window mutations in nadR, pykF, and topA: ", nrow(three.genes) / nrow(filtered.table), "\n")
  
  ## Calc p-value for difference between nadR and pykF+topA in different ways
  
  # t-tests assume normal distribution
  # ks tests look for any difference in the overall distribution (not necessarily an increase)
  
  nadR.window.mutations = filtered.table %>% filter(gene=="nadR")
  pykF.window.mutations = filtered.table %>% filter(gene=="pykF")
  topA.window.mutations = filtered.table %>% filter(gene=="topA")
  
  t.test(nadR.window.mutations$fitness.effect, c(pykF.window.mutations$fitness.effect, topA.window.mutations$fitness.effect) )
  ks.test(nadR.window.mutations$fitness.effect, c(pykF.window.mutations$fitness.effect, topA.window.mutations$fitness.effect) )
  
  t.test(nadR.window.mutations$fitness.effect, topA.window.mutations$fitness.effect)
  ks.test(nadR.window.mutations$fitness.effect, topA.window.mutations$fitness.effect)
  
  t.test(topA.window.mutations$fitness.effect, pykF.window.mutations$fitness.effect)
  ks.test(topA.window.mutations$fitness.effect, pykF.window.mutations$fitness.effect)
  
  # ===> Use these: Tests of medians
  median(nadR.window.mutations$fitness.effect)
  median(topA.window.mutations$fitness.effect)
  median(pykF.window.mutations$fitness.effect)
  wilcox.test(nadR.window.mutations$fitness.effect, topA.window.mutations$fitness.effect, alternative="greater")
  wilcox.test(topA.window.mutations$fitness.effect, pykF.window.mutations$fitness.effect, alternative="greater")
  
  cat("Average fitness effect difference between nadR and topA mutations: ", mean(nadR.window.mutations$fitness.effect) - mean(topA.window.mutations$fitness.effect), "\n")
  cat("Average fitness effect difference between topA and pykF mutations: ", mean(topA.window.mutations$fitness.effect) - mean(pykF.window.mutations$fitness.effect), "\n")
  
  
  # What about those mutations in other genes... how do they compare?
  other.window.mutations = filtered.table %>% filter( !(gene %in% c("nadR", "pykF", "topA")))
  wilcox.test(three.genes$fitness.effect, other.window.mutations$fitness.effect)
  
  
  #######
  
  three.genes.summarize.table = count(three.genes, type, gene)
  three.genes.summarize.table = three.genes.summarize.table %>% group_by(gene) %>% mutate(p=n/sum(n))
  ggplot(three.genes.summarize.table, aes(y = p, x=gene, fill=type), ) + typeColFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0)) 
  ggsave(file.path(summary_output_path,"window_mutation_types_by_gene_scaled.pdf"), width=7.5, height=5)
  
  
  window.population.gene.summarize.table = count(filtered.table, population, gene)
  window.population.gene.summarize.table = window.population.gene.summarize.table %>% group_by(population) %>% mutate(p=n/sum(n))
  
  
  ggplot(filtered.table, aes(x=gene, fill=gene)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + colFillScale + geom_bar(colour="black") + scale_y_continuous(expand=c(0,0), limits = c(0,100)) + scale_x_discrete(drop=FALSE)
  ggsave(file.path(summary_output_path,"window_mutations_by_gene.pdf"), width=8, height=5)
  
  ggplot(window.population.gene.summarize.table, aes(y = p, x=population, fill=gene)) + colFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0)) 
  ggsave(file.path(summary_output_path,"window_by_gene_by_population_scaled.pdf"))
  
  ###############################################################################
  ## Graphs that include only COMPLETE timecourse information
  ###############################################################################
  
  ggplot(complete.table, aes(population, fill=gene)) + colFillScale + geom_bar(colour="black")
  ggsave(file.path(summary_output_path,"complete_by_gene_by_population.pdf"))
  
  complete.population.gene.summarize.table = count(complete.table, population, gene)
  complete.population.gene.summarize.table = complete.population.gene.summarize.table %>% group_by(population) %>% mutate(p=n/sum(n))
  ggplot(complete.population.gene.summarize.table, aes(y = p, x=population, fill=gene)) + colFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0))
  ggsave(file.path(summary_output_path,"complete_by_gene_by_population_scaled.pdf"))
  
  
  ###############################################################################
  ## Graphs for COMBINED timecourse information
  ###############################################################################
  
  combined.three.genes = subset(combined.table, ((gene == "nadR") | (gene == "pykF") | (gene=="topA")))
  
  ## Print the number of combined mutations in each gene
  print(combined.three.genes %>% group_by(gene) %>% summarize(n=n()))
  
  ggplot( combined.three.genes, aes(x=gene, fill=type)) + typeColFillScale + geom_bar(colour="black") + scale_y_continuous(expand = c(0,0))
  ggsave(file.path(summary_output_path,"combined_mutation_types_by_gene.pdf"), width=7.5, height=5)
  
  combined.three.genes.summarize.table = count(combined.three.genes, type, gene)
  combined.three.genes.summarize.table$type = factor(combined.three.genes.summarize.table$type, levels=rev(c("IS_element", "indel", "big_indel", "nonsense", "nonsynonymous", "synonymous", "intergenic")))
  
  combined.three.genes.summarize.table = combined.three.genes.summarize.table %>% group_by(gene) %>% mutate(p=n/sum(n))
  ggplot(combined.three.genes.summarize.table, aes(y = p, x=gene, fill=type), ) + typeColFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0)) 
  ggsave(file.path(summary_output_path,"combined_mutation_types_by_gene_scaled.pdf"), width=7.5, height=5)
  
  ################################################################################################################################
  # Block for counting each **unique** mutation aka unique mutations (counting one time if it occurred in multiple populations!)
  ## We **don't** use this for stats except for exact parallelism at the genetic level.

  combined.unique.table = combined.table %>% group_by(gene, type, position, full_name) %>% summarize(populations = paste(population, collapse=","), num.populations=n())
  write.csv(combined.unique.table, file=file.path(summary_output_path, "complete_all_mutations_unique_events.csv"))

  cat("Number of mutations found in this many populations...\n")
  table(combined.unique.table$num.populations)
  
  ### We now mark which ones are recurrent in the window time course
  
  combined.unique.table.for.join = combined.unique.table %>% ungroup() %>% select(full_name, populations, num.populations)
  window.recurrent.marked.table = window.table %>% left_join(combined.unique.table.for.join, by="full_name")
  window.recurrent.marked.table.recurrent = window.recurrent.marked.table %>% filter(num.populations>1)
  window.recurrent.marked.table.onetime = window.recurrent.marked.table %>% filter(num.populations==1)
  
  write.csv(window.recurrent.marked.table, file=file.path(summary_output_path, "window_all_mutations_unique_events.csv"))
  
  cat("Number of recurrent mutation fitness measurements: ", nrow(window.recurrent.marked.table.recurrent), "\n")
  cat("Number of one-time mutation fitness measurements: ", nrow(window.recurrent.marked.table.onetime), "\n")
  
  ## How many recurrent mutations do we have fitness measurements for?
  ## Each one has the same fitness value = one measurement so collapse this list
  window.recurrent.marked.table.recurrent.measured = window.recurrent.marked.table.recurrent %>% group_by(full_name, populations, num.populations, fitness.effect) %>% summarize(measured.populations = paste(population, collapse=","), num.measured.populations=n())
  
  window.recurrent.marked.table.onetime.measured = window.recurrent.marked.table.onetime
  
  cat("Number of recurrent mutations with fitness measurements: ", nrow(window.recurrent.marked.table.recurrent.measured), "\n")
  
  ## Do the recurrent ones have a higher fitness?
  wilcox.test(window.recurrent.marked.table.recurrent.measured$fitness.effect, window.recurrent.marked.table.onetime.measured$fitness.effect, alternative="greater")
  
  cat("Average fitness effect difference between recurrent and one-time mutations: ", mean(window.recurrent.marked.table.recurrent.measured$fitness.effect) - mean(window.recurrent.marked.table.onetime.measured$fitness.effect), "\n")
  
  
  ################################################################################################################################
  ##
  ## Graphs of each gene showing selection coefficients and files translating positions to within-gene coords
  
  ## This lets us choose what to use here ==> it's all mutations (unique per population)
  for.gene.combined.table = combined.table
  
  
  # Update to mock fitnesses used for graphing
  for.gene.combined.table$fitness.effect[is.na(for.gene.combined.table$fitness.effect)] = 0.01
  for.gene.combined.table$fitness.effect.CI95L[is.na(for.gene.combined.table$fitness.effect.CI95L)] = 0
  for.gene.combined.table$fitness.effect.CI95U[is.na(for.gene.combined.table$fitness.effect.CI95U)] = 0
  
  for.gene.combined.table$gene = factor(for.gene.combined.table$gene, levels=names(myColors))
  
  ggplot(for.gene.combined.table, aes(population, fill=gene)) + colFillScale + geom_bar(colour="black")
  ggsave(file.path(summary_output_path,"complete_mutations_by_gene_by_population.pdf"), width=8, height=5)
  
  population.gene.summarize.table = count(for.gene.combined.table, population, gene)
  population.gene.summarize.table = population.gene.summarize.table %>% group_by(population) %>% mutate(p=n/sum(n))
  ggplot(population.gene.summarize.table, aes(y = p, x=population, fill=gene) ) + colFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0))
  ggsave(file.path(summary_output_path,"complete_mutations_by_gene_by_population_scaled.pdf"), width=8, height=5)
  
  ggplot(for.gene.combined.table, aes(x=gene, fill=gene)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + colFillScale + geom_bar(colour="black") + scale_y_continuous(expand=c(0,0), limits = c(0,120)) + scale_x_discrete(drop=FALSE)
  ggsave(file.path(summary_output_path,"complete_mutations_by_gene.pdf"), width=8, height=5)
  
  for (i in 1:length(levels(for.gene.combined.table$gene))) {
    the.gene = as.character(levels(for.gene.combined.table$gene)[i])
    cat(the.gene, "\n")
    this.table = subset(for.gene.combined.table, gene==the.gene)
    
    if (nrow(this.table) > 0) {
      
      #transform the coordinates so that +1 is the first position of the gene
      gene_start = gene_starts[[the.gene]]
      gene_end = gene_ends[[the.gene]]
      gene_size = gene_sizes[[the.gene]]
      
      cat("Start/End/Size (in base pairs):", gene_start, "/" , gene_end, "/", gene_size, "\n")
      
      if ( gene_strands[[the.gene]] == -1) {
        this.table$position =  gene_size - this.table$position + 1
        gene_start  =  gene_size - gene_end + 1
        gene_end  =  gene_size - gene_start + 1
      } 
      
      this.table$position = this.table$position - gene_start+1
      
      this.table$found.in.both = revalue(this.table$found.in.both, c("1"="one", "2"="both"))
      
      write.csv(this.table, file=file.path(summary_output_path, paste0("combined_graphed_gene_positions_", the.gene, ".csv")))
      
      gene_end = gene_end - gene_start +1
      gene_start = gene_start - gene_start +1
      
      this.with.fitness.effect = subset(this.table, fitness.effect != 0.01)
      mean.fitness.effect = mean(this.with.fitness.effect$fitness.effect)
      
      p = ggplot(this.table, aes(x=position, y=fitness.effect, color=type, shape=found.in.both)) 
      p + geom_hline(aes(yintercept = mean.fitness.effect), color="red", size=0.8, linetype="solid") + geom_vline(aes(xintercept = 1), color="grey", size=0.8, linetype="solid") + geom_vline(aes(xintercept = gene_end), color="grey", size=0.8, linetype="solid") + geom_point(size=3) + geom_errorbar(aes(ymin=fitness.effect.CI95L, ymax=fitness.effect.CI95U), width=0.5) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size=9)) + coord_cartesian(ylim=c(0,0.25), xlim=c(-500,gene_end+500)) + typeColScale
      ggsave(file.path(summary_output_path,filename=paste(the.gene, ".pdf", sep="")), width=10, height=6)
      
      ggplot(this.table, aes(x=fitness.effect, fill=type)) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + coord_cartesian() + typeColFillScale
      ggsave(file.path(summary_output_path,paste("combined_selection_coefficients_by_type_", the.gene, ".pdf", sep="")), width=8, height=5) 
      
      
      ## Test fitness effects of nonsynonymous mutations versus disruptive mutations (nonsense, IS_element, indel, big_indel)
      this.subset.table = subset(window.table, gene==the.gene)
      
      ## We leave out synonymous mutations when making these comparisons
      nonsynonymous_subset = subset(this.subset.table, type=="nonsynonymous")
      disruptive_subset = subset(this.subset.table, type %in% c("nonsense", "IS_element", "indel", "big_indel"))
      res.two.tailed = data.frame(p.value = NA)
      try( {res.two.tailed = wilcox.test(nonsynonymous_subset$fitness.effect, disruptive_subset$fitness.effect) })
      res.greater = data.frame(p.value = NA)
      try( {res.greater = wilcox.test(nonsynonymous_subset$fitness.effect, disruptive_subset$fitness.effect, alternative="greater") })
      res.less = data.frame(p.value = NA)
      try( {res.less = wilcox.test(nonsynonymous_subset$fitness.effect, disruptive_subset$fitness.effect, alternative="less") })
      
      cat (the.gene, " mean of all mutations in gene: ", mean(this.subset.table$fitness.effect), "\n")
      cat (the.gene, " nonsynonymous number: ", nrow(nonsynonymous_subset), " ", " disruptive number:", nrow(disruptive_subset), " ", " other number:", nrow(this.subset.table) - nrow(nonsynonymous_subset) - nrow(disruptive_subset), "\n")
      
      cat (the.gene, " Two tailed p-value of Mann-Whitney U test on nonsynonymous different from disruptive fitness effects: ", res.two.tailed$p.value, "\n")
      cat (the.gene, " One tailed p-value of Mann-Whitney U test on nonsynonymous greater than disruptive fitness effects: ", res.greater$p.value, "\n")
      cat (the.gene, " One tailed p-value of Mann-Whitney U test on disruptive greater than nonsynonymous fitness effects: ", res.less$p.value, "\n")
      

            cat (the.gene, " nonsynonymous mean: ", mean(nonsynonymous_subset$fitness.effect), " ", " disruptive mean:", mean(disruptive_subset$fitness.effect), "\n")
      cat (the.gene, " nonsynonymous median: ", median(nonsynonymous_subset$fitness.effect), " ", " not distruptive median:", median(disruptive_subset$fitness.effect), "\n")
    }
  }
  
  
  ################################################################################################################################
  ### nadR positional stats that use all combined (window+complete) mutations
  
  # Read back in tables to test specific hypotheses about distributions of where the mutations are 
  # in different protein domains
  
  nadR.mutations = read.csv(file.path(summary_output_path, paste0("combined_graphed_gene_positions_",  "nadR", ".csv")))
  
  cat("total mutations within nadR", nrow(nadR.mutations))
  
  #Some of the deletions are so large that we can't localize them to one domain. Remove them before doing stats on domains.
  nadR.mutations = nadR.mutations %>% filter(type != "big_indel")
  
  nadR.end.HTH = 186
  nadR.end.NMN_adenylyltransferase = 687
  
  # the minus 3 bases removes the stop codon
  nadR.size = gene_ends["nadR"] - gene_starts["nadR"] - 3 + 1
  nadR.first.domain = nadR.mutations %>% filter((position <= nadR.end.HTH) & (position >= 1))
  nadR.second.domain = nadR.mutations %>% filter((position > nadR.end.HTH) & (position <= nadR.end.NMN_adenylyltransferase))
  nadR.third.domain = nadR.mutations %>% filter((position > nadR.end.NMN_adenylyltransferase) & (position <= nadR.size))

  cat("nadR mutations within gene and not big indels ", nrow(nadR.first.domain) + nrow(nadR.second.domain) + nrow(nadR.third.domain))
  
  disruptive.mutation.types = c("IS_element", "nonsense", "indel", "big_indel")
  nadR.first.domain.disruptive.mutation = nadR.first.domain %>% filter(type %in% disruptive.mutation.types)
  nadR.first.domain.nondisruptive.mutation = nadR.first.domain %>% filter(!(type %in% disruptive.mutation.types))
  
  nadR.second.domain.disruptive.mutation = nadR.first.domain %>% filter(type %in% disruptive.mutation.types)
  nadR.second.domain.nondisruptive.mutation = nadR.first.domain %>% filter(!(type %in% disruptive.mutation.types))
  
  nadR.third.domain.disruptive.mutation = nadR.third.domain %>% filter(type %in% disruptive.mutation.types)
  nadR.third.domain.nondisruptive.mutation = nadR.third.domain %>% filter(!(type %in% disruptive.mutation.types))
  
  ## This is the test in the text
  fisher.test( matrix(c(nrow(nadR.first.domain.disruptive.mutation) + 
                          nrow(nadR.second.domain.disruptive.mutation),
                        nrow(nadR.first.domain.nondisruptive.mutation) + 
                          nrow(nadR.second.domain.nondisruptive.mutation),
                        nrow(nadR.third.domain.disruptive.mutation),
                        nrow(nadR.third.domain.nondisruptive.mutation)
  ),
  nrow=2, 
  byrow=F
  ),
  alternative="greater"
  )
  
  ################################################################################################################################
  ### pykF positional stats that use all combined (window+complete) mutations
  
  pykF.mutations = read.csv(file.path(summary_output_path, paste0("combined_graphed_gene_positions_", "pykF", ".csv")))
  
  cat("total mutations within pykF", nrow(pykF.mutations))
  
  #Some of the deletions are so large that we can't localize them to one domain. Remove them before doing stats on domains.
  pykF.mutations = pykF.mutations %>% filter(type != "big_indel")
  
  # the minus 3 bases removes the stop codon
  pykF.size = gene_ends["pykF"] - gene_starts["pykF"] - 3 + 1
  pykF.domain.A = pykF.mutations %>% filter( ((position >= 1) & (position <= 210)) | ((position >= 511) & (position <= 1035)))
  pykF.domain.B = pykF.mutations %>% filter((position >= 211) & (position <= 510))
  pykF.domain.C = pykF.mutations %>% filter((position >= 1054) & (position <= 1410))
  
  cat("pykF mutations within gene and not big indels ", nrow(pykF.domain.A) + nrow(pykF.domain.B) +nrow(pykF.domain.C))
  
  
  pykF.domain.A.ns = pykF.domain.A %>% filter(type=="nonsynonymous")
  pykF.domain.B.ns = pykF.domain.B %>% filter(type=="nonsynonymous")
  pykF.domain.C.ns = pykF.domain.C %>% filter(type=="nonsynonymous")
  
  ns.in.domain.A = nrow(pykF.domain.A.ns)
  ns.in.domain.B.or.C = nrow(pykF.domain.B.ns) + nrow(pykF.domain.C.ns)
  length.basis.A = 210 + (1035-511+1)
  length.basis.B.or.C = 510-211+1 + 1410-1054+1
  
  # This is the test in the text
  print(binom.test(ns.in.domain.A, ns.in.domain.A+ns.in.domain.B.or.C, p = length.basis.A/(length.basis.A+length.basis.B.or.C), alternative="greater"))
  
  ################################################################################################################################
  ### topA positional stats that use all combined (window+complete) mutations
  
  topA.mutations = read.csv(file.path(summary_output_path, paste0("combined_graphed_gene_positions_", "topA", ".csv")))
  
  #Some of the deletions are so large that we can't localize them to one domain. Remove them before doing stats on domains.
  topA.mutations = topA.mutations %>% filter(type != "big_indel")
  
  topA.size = gene_ends["topA"] - gene_starts["topA"] - 3 + 1
  
  #These do not include the whole length of the gene
  topA.D1 = topA.mutations %>% filter(((position >= 1) & (position <= 105)) | ((position >= 247) & (position <= 471)))
  topA.D4 = topA.mutations %>% filter(((position >= 106) & (position <= 246)) | ((position >= 472) & (position <= 639)) | ((position >= 1414) & (position <= 1767)))
  topA.catalytic.core = topA.mutations %>% filter(position <= 1767)
  
  topA.D1.ns = topA.D1 %>% filter(type=="nonsynonymous")
  topA.D4.ns = topA.D4 %>% filter(type=="nonsynonymous")
  topA.catalytic.core.ns = topA.catalytic.core %>% filter(type=="nonsynonymous")
  
  ns.in.D1.or.D4 = nrow(topA.D1.ns) + nrow(topA.D4.ns)
  ns.in.catalytic.core = nrow(topA.catalytic.core.ns)
  length.basis.D1.or.D4 = 639-1+1 + (1767-1414+1)
  length.basis.catalytic.core= 1767-1+1
  
  binom.test(ns.in.D1.or.D4, ns.in.catalytic.core, p = length.basis.D1.or.D4/(length.basis.catalytic.core), alternative="greater")
  
}

