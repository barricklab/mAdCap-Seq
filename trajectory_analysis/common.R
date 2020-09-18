library(ggplot2)
library(plyr)
library(grid)
library(dplyr)

autocorrelation_coefficient <- function(v) {
  
  v.mean = mean(v)
  v.ss = 0
  for(i in 1:length(v)) {
    v.ss = v.ss + (v[i]-v.mean)^2
  }
  
  v.num = 0
  for(i in 2:length(v)) {
    v.num = v.num + (v[i]-v.mean) *  (v[i-1]-v.mean)
  }
  
  return(v.num/v.ss)
}


#Evidence mode controls whether RA or JC versus mutation categories are assigned
load_count_files <- function(variant_read_counts_file_name, total_read_counts_file_name, evidence.mode = F) {
  
  # read in the tables, get generations from first columns then nuke that column
  total_read_counts = read.table(total_read_counts_file_name, header=T)
  generations = total_read_counts[,1]
  total_read_counts = total_read_counts[-1]
  variant_read_counts = read.table(variant_read_counts_file_name, header=T)
  variant_read_counts = variant_read_counts[-1]
  
  ## Set up a main table with information about each mutation
  full_mutation_names = colnames(total_read_counts)
  split_names = strsplit(full_mutation_names, ".", fixed=T)
  i=1
  length(split_names)
  
  mutation_info = data.frame()
  for(j in 1:length(split_names))  {
    line = split_names[[j]]
    
    if (line[2] == "JC") {
      
      this.category = "JC";
      this.type = "junction"
      
      this.position = as.numeric(line[6])
      this.start.position = this.position
      this.end.position = this.start.position
      
      #this indicates a rather short indel (not IS-mediated)
      if (line[1] == line[5]) {
        this.type = "indel"
        
        this.start.position = min(as.numeric(line[3]), as.numeric(line[6]))
        this.end.position = max(as.numeric(line[3]), as.numeric(line[6]))
        
        if (this.end.position - this.start.position > 50) {
          this.type = "big_indel"
        }
      }
            
      new_mutation_info = data.frame(gene=line[5], full_name=full_mutation_names[[i]], type=this.type, category=this.category, position=this.position, start.position=this.start.position, end.position=this.end.position, freq.corr=as.numeric(line[9]))
      mutation_info = rbind(mutation_info, new_mutation_info)  
      
    } else if (line[2] == "RA") {
      
      this.category = "RA"
      
      this.type = "unknown"
      if (line[9] == "snp_type") {
        this.type = line[10]
      } else if (line[10] == "snp_type") {
        this.type = line[11]
      }
      
      #get rid of 'intergenic' assignments overriding indels
      if ((line[4] == "_") || (line[5] == "_")) {
        this.type = "indel"
      }
      
      this.position = as.numeric(line[3])
            
      new_mutation_info = data.frame(gene=line[1], full_name=full_mutation_names[[i]], type=this.type, category=this.category, position=this.position, start.position=this.position, end.position=this.position, freq.corr=1)
      mutation_info = rbind(mutation_info, new_mutation_info) 
      
    } else if (line[2] == "MOB") { 
      
      this.position = as.numeric(line[6])
      this.duplication.size = as.numeric(line[7])
      this.start.position = this.position
      this.end.position = this.start.position + this.duplication.size - 1
      
      new_mutation_info = data.frame(gene=line[5], full_name=full_mutation_names[[i]], type="IS_element", category="MOB", position=this.position, start.position=this.start.position, end.position=this.end.position, freq.corr = 1 )
      mutation_info = rbind(mutation_info, new_mutation_info) 
      
    } else {
      # Unknown = delete!
      stop(paste("unknown line type", line[2]))
      total_read_counts = total_read_counts[-i]
      variant_read_counts = variant_read_counts[-i]
      i = i-1
    }
    i = i+1
  }
  
  #Fix the names
  mutation_info$gene = revalue(mutation_info$gene, c("REL606_hslu"="hslU", "REL606_iclr"="iclR", "REL606_nadr"="nadR", "REL606_pykf"="pykF", "REL606_spot"="spoT", "REL606_topa"="topA", "REL606_ybal"="ybaL", "REL606_yijc"="fabR", "REL606"="other"))
  
  ## merge to create an unrolled table for fitting and graphing
  big_num_rows = nrow(mutation_info)*nrow(variant_read_counts)
  big.full_names = character(big_num_rows)
  big.genes = character(big_num_rows)
  big.types = character(big_num_rows)
  big.evidences = numeric(big_num_rows)
  big.generations = numeric(big_num_rows)
  big.variant_read_counts = numeric(big_num_rows)
  big.total_read_counts = numeric(big_num_rows)
  big.L95.CIs = numeric(big_num_rows)
  big.U95.CIs = numeric(big_num_rows)
  big.freq.corrs = numeric(big_num_rows)
  
  
  k = 0
  for(j in 1:nrow(mutation_info))  {
    for (i in 1:nrow(variant_read_counts)) {
      k = k + 1
      big.full_names[k] = as.character(mutation_info$full_name[j])
      big.genes[k] = as.character(mutation_info$gene[j])
      big.types[k] = as.character(mutation_info$type[j])
      big.evidences[k]=j
      big.generations[k] = generations[i]
      big.variant_read_counts[k] = variant_read_counts[i,j]
      big.total_read_counts[k] = total_read_counts[i,j]
      big.L95.CIs[k] = NA
      big.U95.CIs[k] = NA
      big.freq.corrs[k] = mutation_info$freq.corr[j]
    }
  }
  
  big_table = data.frame(full_name = as.factor(big.full_names), gene = as.factor(big.genes), type = as.factor(big.types), evidence = big.evidences, generation = big.generations, variant_read_count = big.variant_read_counts, total_read_count = big.total_read_counts, L95.CI = big.L95.CIs, U95.CI = big.U95.CIs, freq.corr = big.freq.corrs)
  
  #Calculate some additional columns
  big_table$frequency = (big_table$variant_read_count / big_table$total_read_count) * big_table$freq.corr
  big_table$log_frequency = log(big_table$frequency, base=10)
  big_table$not_variant_read_count = big_table$total_read_count - big_table$variant_read_count
  
  return(list(mutation_info=mutation_info, big_table=big_table))
}