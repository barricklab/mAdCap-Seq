combination_graphs <- function () 
{
  
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(scales) 
  
  output_path = "grouped_output"
  summary_output_path = file.path(output_path, "summary")
  dir.create(summary_output_path, showWarnings = FALSE)
  
  
  file_names=c("A1_population_window_consensus_straightened_significant_mutations.csv", "A2_population_window_consensus_straightened_significant_mutations.csv", "A3_population_window_consensus_straightened_significant_mutations.csv", "A6_population_window_consensus_straightened_significant_mutations.csv", "A7_population_window_consensus_straightened_significant_mutations.csv", "A9_population_window_consensus_straightened_significant_mutations.csv")
  
  #file_names=c("A1_population_window_straightened_significant_mutations.csv", "A2_population_window_straightened_significant_mutations.csv", "A3_population_window_straightened_significant_mutations.csv", "A6_population_window_straightened_significant_mutations.csv", "A7_population_window_straightened_significant_mutations.csv", "A9_population_window_straightened_significant_mutations.csv")
  
  names=c("A1", "A2", "A3", "A6", "A7", "A9")
  
  myTypeColors <- c("orange", "black", "gray", "blue", "red", "green", "purple", "magenta")
  names(myTypeColors) <- c("IS_element", "indel", "big_indel", "nonsynonymous", "nonsense", "synonymous","intergenic", "junction")
  typeColScale <- scale_colour_manual(name = "type",values = myTypeColors)
  typeColFillScale = scale_fill_manual(name = "gene",values = myTypeColors)
  
  
  gene_sizes = c( pykF=4443, fabR=3781, ybaL=4533, hslU=4738, nadR=4032, iclR=3710, spoT=5231, topA=5656)
  gene_strands = c( pykF=+1, fabR=+1, ybaL=-1, hslU=-1, nadR=+1, iclR=-1, spoT=+1, topA=5656)
  gene_starts = c( pykF=1631, fabR=1677, ybaL=1400, hslU=1400, nadR=1400, iclR=1400, spoT=1429, topA=1659)
  gene_ends = c( pykF=3043, fabR=2381, ybaL=3076, hslU=2731, nadR=2632, iclR=2224, spoT=3831, topA=4256)
  
  #Setup for graphs
  theme_set(theme_bw(base_size = 24))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
  myColors <- c("green", "orange", "blue", "red", "cyan", "brown", "purple", "magenta", "grey")
  names(myColors) <- c("nadR", "pykF", "topA", "spoT", "fabR", "ybaL", "hslU", "iclR", "other")
  colScale <- scale_colour_manual(name = "Gene",values = myColors)
  colFillScale = scale_fill_manual(name = "Gene",values = myColors)
  
  combined.table = data.frame()
  for (i in 1:length(file_names)) {
    
    this.table = read.csv(file.path(output_path, "output", file_names[i]), header=T)
    this.table$population = names[i]
    combined.table = rbind(combined.table, this.table)
  }
  combined.table$full_name_population = paste(combined.table$full_name,combined.table$population , sep=".")
  combined.table$gene = factor(as.character(combined.table$gene), levels=names(myColors))
  combined.table$type = factor(as.character(combined.table$type), levels=names(myTypeColors))
  
  write.csv(combined.table, file.path(output_path, "output", paste("consensus_significant_mutations.csv", sep="")), row.names=F)
  
  nrow(combined.table)
  filtered.table = combined.table
  #filtered.table = subset(combined.table, fitness.effect<0.09)
  nrow(filtered.table)
  
  this.max.y = 25
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=gene)) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + colFillScale + scale_y_continuous(expand = c(0,0), limits = c(0,this.max.y))
  ggsave(file.path(summary_output_path,"fitness_effects_by_gene.pdf"), width=8, height=5)
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=gene)) + geom_histogram(binwidth=0.0025, alpha=1, position="stack", boundary = 0.5) + colFillScale + scale_y_continuous(expand = c(0,0), limits = c(0,16)) + facet_grid(rows=vars(gene))
  ggsave(file.path(summary_output_path,"fitness_effects_by_gene.pdf"), width=8, height=5)
  
  
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=type)) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + scale_y_continuous(expand = c(0,0), limits = c(0,this.max.y)) + typeColFillScale
  ggsave(file.path(summary_output_path,"fitness_effects_by_type.pdf"), width=8, height=5) 
  
  facet_table = subset(filtered.table, type %in% c("IS_element", "indel", "nonsynonymous"))
  facet_table = subset(facet_table, gene %in% c("nadR", "pykF", "topA"))
  ggplot(facet_table, aes(x=fitness.effect, fill=population)) +facet_grid(type~gene) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + scale_y_continuous(expand = c(0,0)) 
  ggsave(file.path(summary_output_path,"selection_coefficients_by_type_by_gene.pdf"), width=8, height=5) 
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=population)) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + scale_y_continuous(expand = c(0,0))
  ggsave(file.path(summary_output_path,"selection_coefficients_by_population.pdf"), width=8, height=5)
  
  ggplot(filtered.table, aes(x=fitness.effect, fill=type)) + geom_histogram(binwidth=0.002, alpha=1, position="stack", boundary=0) + scale_y_continuous(expand = c(0,0)) + facet_wrap(~population) + typeColFillScale
  ggsave(file.path(summary_output_path,"selection_coefficients_by_population_by_type_colors.pdf"), width=8, height=5)
  
  ggplot(filtered.table, aes(population, fill=gene)) + colFillScale + geom_bar(colour="black") + scale_y_continuous(expand = c(0,0))
  ggsave(file.path(summary_output_path,"selection_coefficients_fit_by_gene_by_population.pdf"), width=8, height=5)
  
  
  ####### stats interlude
  
  ## Calc average and sd of all selection coefficients
  mean(filtered.table$fitness.effect)
  sd(filtered.table$fitness.effect)
  
  
  #######
  # How many times did we see the same mutation in different populations?
  
  mutations.in.populations = filtered.table %>% group_by(full_name) %>% summarize(num.populations=n(), populations=paste0(population, collapse = "-")) 
  
  mutations.in.multiple.populations = mutations.in.populations %>% filter(num.populations > 1)
  
  write.csv(mutations.in.multiple.populations, file.path(summary_output_path,"mutations_in_multiple_populations.csv"))
  
  cat("Number of mutations that were found in more than one population: ", nrow(mutations.in.multiple.populations))
  
  cat("Number of mutations (counting multi-populations multiple times): ", nrow(filtered.table))
  cat("Number of distinct mutations (counting multi-populations once): ", nrow(filtered.table) - sum(mutations.in.multiple.populations$num.populations-1 ))
  
  
  ## Make the table look nice here....
  filtered.table %>% left_join(mutations.in.populations, by="full_name")
  
  write.csv(filtered.table, file.path(output_path, "output", paste("consensus_significant_mutations.csv", sep="")), row.names=F)
  
  
  #############
  # Is there a difference in the distribution of fitness effects across populations?
  
  kruskal.test(fitness.effect ~ population, data=filtered.table)
  
  # Yes, but it could be for stochastic reasons on the high side and
  # sampling depth on the low side
  
  ggplot(filtered.table, aes(x=population, y=fitness.effect)) + 
    #colFillScale + 
    geom_violin(trim=FALSE) + 
    geom_jitter(width=0.2)
  
  ggsave(file.path(summary_output_path,"violin_fitness_effects_by_population.pdf"), width=8, height=5)
  
  ##############
  
  three.genes = subset(filtered.table, ((gene == "nadR") | (gene == "pykF") | (gene=="topA")))
  ggplot( three.genes, aes(x=gene, fill=type)) + typeColFillScale + geom_bar(colour="black") + scale_y_continuous(expand = c(0,0))
  ggsave(file.path(summary_output_path,"mutation_types_by_gene.pdf"), width=7.5, height=5)
  
  
  ## Calc p-value for difference between nadR and pykF+topA
  nadR.window.mutations = filtered.table %>% filter(gene=="nadR")
  pykF.window.mutations = filtered.table %>% filter(gene=="pykF")
  topA.window.mutations = filtered.table %>% filter(gene=="topA")
  
  t.test(nadR.window.mutations$fitness.effect, c(pykF.window.mutations$fitness.effect, topA.window.mutations$fitness.effect) )
  ks.test(nadR.window.mutations$fitness.effect, c(pykF.window.mutations$fitness.effect, topA.window.mutations$fitness.effect) )
  
  t.test(nadR.window.mutations$fitness.effect, topA.window.mutations$fitness.effect)
  ks.test(nadR.window.mutations$fitness.effect, topA.window.mutations$fitness.effect)
  
  t.test(topA.window.mutations$fitness.effect, pykF.window.mutations$fitness.effect)
  ks.test(topA.window.mutations$fitness.effect, pykF.window.mutations$fitness.effect)
  
  #alternative tests of medians
  median(nadR.window.mutations$fitness.effect)
  median(topA.window.mutations$fitness.effect)
  median(pykF.window.mutations$fitness.effect)
  wilcox.test(nadR.window.mutations$fitness.effect, topA.window.mutations$fitness.effect, alternative="greater")
  wilcox.test(topA.window.mutations$fitness.effect, pykF.window.mutations$fitness.effect, alternative="greater")
  
  
  #######
  
  three.genes.summarize.table = count(three.genes, type, gene)
  three.genes.summarize.table = three.genes.summarize.table %>% group_by(gene) %>% mutate(p=n/sum(n))
  ggplot(three.genes.summarize.table, aes(y = p, x=gene, fill=type), ) + typeColFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0)) 
  ggsave(file.path(summary_output_path,"mutation_types_by_gene_scaled.pdf"), width=7.5, height=5)
  
  
  window.population.gene.summarize.table = count(filtered.table, population, gene)
  window.population.gene.summarize.table = window.population.gene.summarize.table %>% group_by(population) %>% mutate(p=n/sum(n))
  
  
  ggplot(filtered.table, aes(x=gene, fill=gene)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + colFillScale + geom_bar(colour="black") + scale_y_continuous(expand=c(0,0), limits = c(0,100)) + scale_x_discrete(drop=FALSE)
  ggsave(file.path(summary_output_path,"window_mutations_by_gene.pdf"), width=8, height=5)
  
  ggplot(window.population.gene.summarize.table, aes(y = p, x=population, fill=gene)) + colFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0)) 
  ggsave(file.path(summary_output_path,"window_by_gene_by_population_scaled.pdf"))
  
  ## full population stuff
  
  file_names=c("A1_population_complete_significant_mutations.csv", "A2_population_complete_significant_mutations.csv", "A3_population_complete_significant_mutations.csv", "A7_population_complete_significant_mutations.csv")
  names=c("A1", "A2", "A3", "A7")
  
  full.combined.table = data.frame()
  for (i in 1:length(file_names)) {
    
    this.table = read.csv(file.path(output_path, "output", file_names[i]), header=T)
    this.table$population = names[i]
    full.combined.table = rbind(full.combined.table, this.table)
  }
  full.combined.table$full_name_population = paste(full.combined.table$full_name,full.combined.table$population , sep=".")
  
  ggplot(full.combined.table, aes(population, fill=gene)) + colFillScale + geom_bar(colour="black")
  ggsave(file.path(summary_output_path,"complete_by_gene_by_population.pdf"))
  
  complete.population.gene.summarize.table = count(full.combined.table, population, gene)
  complete.population.gene.summarize.table = complete.population.gene.summarize.table %>% group_by(population) %>% mutate(p=n/sum(n))
  ggplot(complete.population.gene.summarize.table, aes(y = p, x=population, fill=gene)) + colFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0))
  ggsave(file.path(summary_output_path,"complete_by_gene_by_population_scaled.pdf"))
  
  union_table = table(c(as.character(full.combined.table$full_name_population), as.character(filtered.table$full_name_population)))
  
  
  full.combined.table$fitness.effect = 0.01
  full.combined.table$fitness.effect.CI95L = 0.0
  full.combined.table$fitness.effect.CI95U = 0.0
  
  #full.combined.table = data.frame()
  
  full.combined.table$found.in.both = as.numeric(union_table[full.combined.table$full_name_population])
  
  
  
  
  full.combined.table = subset (full.combined.table, found.in.both==1)
  
  
  for.gene.combined.table = data.frame(gene = as.factor(c(as.character(full.combined.table$gene), as.character(filtered.table$gene))), fitness.effect = c(full.combined.table$fitness.effect, filtered.table$fitness.effect), fitness.effect.CI95L = c(full.combined.table$fitness.effect.CI95L, filtered.table$fitness.effect.CI95L), fitness.effect.CI95U = c(full.combined.table$fitness.effect.CI95U, filtered.table$fitness.effect.CI95U), type = as.factor(c(as.character(full.combined.table$type), as.character(filtered.table$type))), population = c(full.combined.table$population, filtered.table$population), position = c(full.combined.table$position, filtered.table$position), full_name_population = as.factor(c(as.character(full.combined.table$full_name_population), as.character(filtered.table$full_name_population))) )
  
  for.gene.combined.table$found.in.both  = as.factor(union_table[for.gene.combined.table$full_name_population])
  
  for.gene.combined.table$gene = factor(as.character(for.gene.combined.table$gene), levels=names(myColors[1:8]))
  
  ggplot(for.gene.combined.table, aes(population, fill=gene)) + colFillScale + geom_bar(colour="black")
  ggsave(file.path(summary_output_path,"all_mutations_by_gene_by_population.pdf"), width=8, height=5)
  
  population.gene.summarize.table = count(for.gene.combined.table, population, gene)
  population.gene.summarize.table = population.gene.summarize.table %>% group_by(population) %>% mutate(p=n/sum(n))
  ggplot(population.gene.summarize.table, aes(y = p, x=population, fill=gene) ) + colFillScale + geom_bar(colour="black", stat="identity") + scale_y_continuous(labels = percent_format(), expand = c(0,0))
  ggsave(file.path(summary_output_path,"all_mutations_by_gene_by_population_scaled.pdf"), width=8, height=5)
  
  ggplot(for.gene.combined.table, aes(x=gene, fill=gene)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + colFillScale + geom_bar(colour="black") + scale_y_continuous(expand=c(0,0), limits = c(0,120)) + scale_x_discrete(drop=FALSE)
  ggsave(file.path(summary_output_path,"all_mutations_by_gene.pdf"), width=8, height=5)
  
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
      
      write.csv(this.table, file=file.path(summary_output_path, paste0(the.gene, "_summary_graphed_gene_positions.csv")))
      
      
      gene_end = gene_end - gene_start +1
      gene_start = gene_start - gene_start +1
      
      this.with.fitness.effect = subset(this.table, fitness.effect != 0.01)
      mean.fitness.effect = mean(this.with.fitness.effect$fitness.effect)
      
      p = ggplot(this.table, aes(x=position, y=fitness.effect, color=type, shape=found.in.both)) 
      p + geom_segment(aes(x = -500, y = mean.fitness.effect, xend = gene_end+500, yend = mean.fitness.effect), color="red", size=0.8, linetype="solid") + geom_segment(aes(x = 1, y = 0, xend = 1, yend = 0.1), color="grey", size=0.8, linetype="solid") + geom_segment(aes(x = gene_end, y = 0, xend = gene_end, yend = 0.1), color="grey", size=0.8, linetype="solid") + geom_point(size=3) + geom_errorbar(aes(ymin=fitness.effect.CI95L, ymax=fitness.effect.CI95U), width=0.5) + theme(legend.position="bottom", legend.title=element_blank(), legend.text = element_text(size=9)) + coord_cartesian(ylim=c(0,0.25), xlim=c(-500,gene_end+500)) + typeColScale
      ggsave(file.path(summary_output_path,filename=paste(the.gene, ".pdf", sep="")), width=10, height=6)
      
      ggplot(this.table, aes(x=fitness.effect, fill=type)) + geom_histogram(binwidth=0.002, alpha=1, position="stack") + coord_cartesian() + typeColFillScale
      ggsave(file.path(summary_output_path,paste(the.gene, "_selection_coefficients_by_type.pdf", sep="")), width=8, height=5) 
      
      
      ##nonsynonymous versus other test 
      this.subset.table = subset(this.table, ! ((found.in.both=="one") &(fitness.effect==0.01)))
      
      nonsynonymous_subset = subset(this.subset.table, type=="nonsynonymous")
      not_nonsynonymous_subset = subset(this.subset.table, type!="nonsynonymous")
      res = data.frame(p.value=NA)
      try( {res = wilcox.test(nonsynonymous_subset$fitness.effect, not_nonsynonymous_subset$fitness.effect) })
      
      cat (the.gene, " p-value of Mann-Whitney U test on nonsynonymous different selection coefficients: ", res$p.value, "\n")
      cat (the.gene, " nonsynonymous median: ", median(nonsynonymous_subset$fitness.effect), " ", " not nonsynonymous median:", median(not_nonsynonymous_subset$fitness.effect), "\n")
    }
  }
  
  write.csv(for.gene.combined.table, file=file.path(summary_output_path, "summary_all_mutations.csv"))
  
  ##How many total in the three main genes?
  nrow(for.gene.combined.table %>% filter(gene %in% c("nadR", "pykF", "topA")))
  
  
  ## Read back in tables to test specific hypotheses about distributions in different domains
  
  nadR.mutations = read.csv(file.path(summary_output_path, paste0("nadR", "_summary_graphed_gene_positions.csv")))
  
  nadR.middle = (gene_ends[["nadR"]] - gene_starts[["nadR"]])/2
  nadR.end.HTH = 186
  nadR.end.NMN_adenylyltransferase = 687
  
  nadR.first.domain = nadR.mutations %>% filter(position <= nadR.end.HTH)
  nadR.second.domain = nadR.mutations %>% filter((position > nadR.end.HTH) & (position <= nadR.end.NMN_adenylyltransferase))
  nadR.third.domain = nadR.mutations %>% filter(position > nadR.end.NMN_adenylyltransferase)
  
  nrow(nadR.first.domain) + nrow(nadR.second.domain) + nrow(nadR.third.domain)
  
  disruptive.mutation.types = c("IS_element", "nonsense", "indel", "big_indel")
  nadR.first.domain.disruptive.mutation = nadR.first.domain %>% filter(type %in% disruptive.mutation.types)
  nadR.first.domain.nondisruptive.mutation = nadR.first.domain %>% filter(!(type %in% disruptive.mutation.types))
  
  nadR.second.domain.disruptive.mutation = nadR.first.domain %>% filter(type %in% disruptive.mutation.types)
  nadR.second.domain.nondisruptive.mutation = nadR.first.domain %>% filter(!(type %in% disruptive.mutation.types))
  
  nadR.third.domain.disruptive.mutation = nadR.third.domain %>% filter(type %in% disruptive.mutation.types)
  nadR.third.domain.nondisruptive.mutation = nadR.third.domain %>% filter(!(type %in% disruptive.mutation.types))
  
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
  
  fisher.test( matrix(c(nrow(nadR.first.domain.disruptive.mutation),
                        nrow(nadR.first.domain.nondisruptive.mutation),
                        nrow(nadR.third.domain.disruptive.mutation) + 
                          nrow(nadR.second.domain.disruptive.mutation),
                        nrow(nadR.third.domain.nondisruptive.mutation) + 
                          nrow(nadR.second.domain.nondisruptive.mutation)
  ),
  nrow=2, 
  byrow=F
  )
  )
  
  ks.test(c(nadR.first.domain.disruptive.mutation$fitness.effect,
            nadR.second.domain.disruptive.mutation$fitness.effect,
            nadR.third.domain.disruptive.mutation$fitness.effect
  ),
  c(nadR.first.domain.nondisruptive.mutation$fitness.effect,
    nadR.second.domain.nondisruptive.mutation$fitness.effect,
    nadR.third.domain.nondisruptive.mutation$fitness.effect
  ),
  alternative="greater"
  )
  
  
  
  #######
  
  pykF.mutations = read.csv(file.path(summary_output_path, paste0("pykF", "_summary_graphed_gene_positions.csv")))
  
  pykF.domain.A = pykF.mutations %>% filter((position <= 210) | ((position >= 511) & (position <= 1035)))
  pykF.domain.B = pykF.mutations %>% filter((position >= 211) & (position <= 510))
  pykF.domain.C = pykF.mutations %>% filter((position >= 1054) & (position <= 1410))
  
  nrow(pykF.domain.A) + nrow(pykF.domain.B) +nrow(pykF.domain.C)
  
  pykF.domain.A.ns = pykF.domain.A %>% filter(type=="nonsynonymous")
  pykF.domain.B.ns = pykF.domain.B %>% filter(type=="nonsynonymous")
  pykF.domain.C.ns = pykF.domain.C %>% filter(type=="nonsynonymous")
  
  ns.in.domain.A = nrow(pykF.domain.A.ns)
  ns.in.domain.B.or.C = nrow(pykF.domain.B.ns) + nrow(pykF.domain.C.ns)
  length.basis.A = 210 + (1035-511+1)
  length.basis.B.or.C = 510-211+1 + 1410-1054+1
  
  print(binom.test(ns.in.domain.A, ns.in.domain.A+ns.in.domain.B.or.C, p = length.basis.A/(length.basis.A+length.basis.B.or.C), alternative="greater"))
  
  
  #################
  
  topA.mutations = read.csv(file.path(summary_output_path, paste0("topA", "_summary_graphed_gene_positions.csv")))
  
  topA.D1 = topA.mutations %>% filter(((position >= 1) & (position <= 105)) | ((position >= 247) & (position <= 471)))
  topA.D4 = topA.mutations %>% filter(((position >= 106) & (position <= 246)) | ((position >= 472) & (position <= 639)) | ((position >= 1414) & (position <= 1767)))
  topA.catalytic.core = topA.mutations %>% filter(position <= 1767)
  
  nrow(topA.D1) + nrow(topA.D4)
  
  topA.D1.ns = topA.D1 %>% filter(type=="nonsynonymous")
  topA.D4.ns = topA.D4 %>% filter(type=="nonsynonymous")
  topA.catalytic.core.ns = topA.catalytic.core %>% filter(type=="nonsynonymous")
  
  ns.in.D1.or.D4 = nrow(topA.D1.ns) + nrow(topA.D4.ns)
  ns.in.catalytic.core = nrow(topA.catalytic.core)
  length.basis.D1.or.D4 = 639-1+1 + (1767-1414+1)
  length.basis.catalytic.core= 1767-1+1
  
  binom.test(ns.in.D1.or.D4, ns.in.catalytic.core, p = length.basis.D1.or.D4/(length.basis.catalytic.core), alternative="greater")
  
}
