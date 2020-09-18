grouped_mutations_find <- function (variant_read_counts_file_name, total_read_counts_file_name, output_path, base_file_name)
{
  cat(variant_read_counts_file_name, " ", total_read_counts_file_name, " ", output_path, " ", base_file_name, "\n")
  
  source("common.R")
  library(pvclust)
  
  ### Global settings
  max_distance_betweens = c()
  max_distance_betweens[["RA"]] = 6
  max_distance_betweens[["JC"]] = 20
  
  similarity_thresholds = c()
  similarity_thresholds[["RA"]] = 0.1;
  similarity_thresholds[["JC"]] = 0.5;
  minimum_mutation_frequency = 0.0003
  
  
  exclude_generations = c()
  if (base_file_name=="A1_population_window") {
    exclude_generations = c(163.4, 190.08, 216.76, 243.44)
  } else if (base_file_name=="A2_population_window") {
    exclude_generations = c(163.4, 190.08, 216.76, 243.44)
  } else if (base_file_name=="A3_population_window") {
    exclude_generations = c(163.4, 190.08, 216.76, 243.44)
  } else if (base_file_name=="A6_population_window") {
    exclude_generations = c(163.4, 190.08, 216.76, 243.44)
  } else if (base_file_name=="A7_population_window") {
    exclude_generations = c(163.4, 190.08, 216.76, 243.44)
  } else if (base_file_name=="A9_population_window") {
    exclude_generations = c(163.4, 190.08, 216.76, 243.44)
  } else if (base_file_name=="A1_population_complete") {
    exclude_generations = c(169,183,196,209,223,236)
  } else if (base_file_name=="A2_population_complete") {
    exclude_generations = c(169,183,196,209,223,236,270.12,510.24)
  } else  if (base_file_name=="A3_population_complete") {
    exclude_generations = c(169,183,196,209,223,236)
  } else  if (base_file_name=="A7_population_complete") {
    exclude_generations = c(169,183,196,209,223,236)
  }
  
  #Setup for graphs
  theme_set(theme_bw(base_size = 24))
  theme_update(panel.border=element_rect(color="black", fill=NA, size=1), legend.key=element_rect(color=NA, fill=NA))
  line_thickness = 0.8
  log_frequency_bounds = c(-5,0)
  myColors <- c("purple", "magenta", "orange", "green", "red", "blue", "brown", "cyan", "grey")
  names(myColors) <- c("hslU", "iclR", "pykF", "nadR", "spoT", "topA", "ybaL", "fabR", "other")
  colScale <- scale_colour_manual(name = "Gene",values = myColors)
  
  
  dir.create(file.path(output_path, "graphs"), showWarnings = FALSE)
  dir.create(file.path(output_path, "output"), showWarnings = FALSE)
  
  output_table = data.frame()
  group_counter = 1
  
  res = load_count_files(variant_read_counts_file_name, total_read_counts_file_name)
  
  #total read counts are read in as matrices
  total_read_counts = read.table(total_read_counts_file_name, header=T)
  generations = total_read_counts[,1]
  total_read_counts = total_read_counts[-1]
  variant_read_counts = read.table(variant_read_counts_file_name, header=T)
  variant_read_counts = variant_read_counts[-1]
  
  total_read_counts = total_read_counts[!(generations %in% exclude_generations),]
  variant_read_counts = variant_read_counts[!(generations %in% exclude_generations),]
  
  frequency_matrix = t(as.matrix(variant_read_counts) / as.matrix(total_read_counts))
  rownames(frequency_matrix) = colnames(total_read_counts)
  
  #only care about mutations that reach above some threshold frequency
  
  #Use Canberra distance
  # http://gedas.bizhat.com/dist.htm
  
  for (evidence_type in c("RA","JC")) {
    
    
    
    max_distance_between_linked_mutations = max_distance_betweens[[evidence_type]]
    similarity_threshold = similarity_thresholds[[evidence_type]]
    
    mutation_info <- res$mutation_info 
    mutation_info$min_mutation_frequency = apply(frequency_matrix, 1, max)
    mutation_info = subset(mutation_info, category==evidence_type)
    nrow(mutation_info)
    mutation_info= subset(mutation_info, min_mutation_frequency >= minimum_mutation_frequency)
    nrow(mutation_info)
    
    gene_list = unique(mutation_info$gene)
    
    for (on_gene in gene_list) {
      
      
      mutations_in_gene = subset(mutation_info, gene==on_gene)
      
      
      #sort by positions within gene
      mutations_in_gene = mutations_in_gene[order(mutations_in_gene$position) , ]
      
      #determine which mutations are close enough to compare
      first_mutation_index = 1
      first_position = mutations_in_gene$position[first_mutation_index]
      while (first_mutation_index < nrow(mutations_in_gene) ) {
        last_mutation_index = first_mutation_index
        while ( (last_mutation_index+1 < nrow(mutations_in_gene)) && (mutations_in_gene$position[last_mutation_index+1] < first_position + max_distance_between_linked_mutations)) {
          last_mutation_index = last_mutation_index + 1
        }
        cat( first_mutation_index, " ", last_mutation_index, "\n")
        
        #calculate a distance matrix for these mutations
        if (first_mutation_index != last_mutation_index) {
          on_full_names = as.vector(mutations_in_gene$full_name[first_mutation_index:last_mutation_index])
          
          cat( on_full_names, "\n")
          
          on_total_read_counts = total_read_counts[,on_full_names]
          on_variant_read_counts = variant_read_counts[,on_full_names]
          
          on_frequency_matrix = t(as.matrix(on_variant_read_counts) / as.matrix(on_total_read_counts))
          rownames(on_frequency_matrix) = on_full_names
          
          on_row = 1
          while (on_row <= nrow(on_frequency_matrix) ) {
            if (sum(on_frequency_matrix[on_row,], na.rm = TRUE) == 0) {
              on_frequency_matrix = on_frequency_matrix[-on_row,,drop=FALSE]
            } else {
              on_row = on_row + 1
            }
          }      
          #on_nonzero_frequency_matrix=on_frequency_matrix[!!rowSums(on_frequency_matrix, na.rm = TRUE),drop=FALSE]
          
          #cat("Nonzero rows = ", nrow(on_frequency_matrix) , "\n")
          
          
          if (!is.null(nrow(on_frequency_matrix)) && (nrow(on_frequency_matrix)>1) ) {
            
            cat("Nonzero rows = ", nrow(on_frequency_matrix) , "\n")
            d = dist(on_frequency_matrix, method="canberra") / ncol(on_frequency_matrix)
            
            hc = hclust(d, method="single")
            myhcl <- cutree(hc, h=similarity_threshold)
            print(myhcl)
            
            for (n in 1:max(myhcl)) {
              group_full_names = c()
              for (p in 1:length(myhcl)) {
                if (myhcl[p] == n) {
                  group_full_names = c(group_full_names, names(myhcl)[p])
                }
              }
              
              if (length(group_full_names) > 1) {
                cat("Groups = ", length(unique(myhcl)) , "\n")
                cat("Min Dist = ", min(d) , "\n")
                
                
                group_total_read_counts = total_read_counts[,group_full_names]
                group_variant_read_counts = variant_read_counts[,group_full_names]
                group_frequency_matrix = t(as.matrix(group_variant_read_counts) / as.matrix(group_total_read_counts))
                d = dist(group_frequency_matrix, method="canberra") / ncol(group_frequency_matrix)
                
                on_graph_file_name = paste("groups_", on_gene, "_", evidence_type, "_", first_mutation_index, "-", last_mutation_index, "_group=", group_counter, "_d=", max(d), ".pdf", sep="")
                
                filtered_table = subset(res$big_table, full_name %in% group_full_names)
                filtered_table = subset(filtered_table, !(generation %in% exclude_generations))
                
                p = ggplot(filtered_table, aes(x=generation, y=frequency, color=full_name, group=evidence))
                p + geom_line(size=line_thickness) + theme(legend.text = element_text(size=5))
                ggsave(filename=file.path(output_path, "graphs", base_file_name, paste(on_graph_file_name, sep=""))) 
                
                
                for (q in 1:length(group_full_names)) {
                  new_row = data.frame(full_name=group_full_names[q], group=group_counter, evidence=evidence_type, max_distance = max(d)) 
                  output_table = rbind(output_table, new_row)
                }
                
                
                group_counter = group_counter + 1
              }
              
            }
          }
        }
        first_mutation_index = last_mutation_index + 1
        first_position = mutations_in_gene$position[first_mutation_index]
      }
    }
  }
  
  on_graph_file_name = file.path(output_path, "output", paste(base_file_name, "_grouped_mutations.csv", sep=""))
  write.csv(output_table, on_graph_file_name)
  
}
