source("common.R")
library(dplyr)
options(dplyr.summarise.inform=F) 
library(forcats)
library(stringr)

LTEE_compare <- function() {

  full.table = read.csv(file.path("LTEE-input/LTEE_mutation_comparison.csv"))

  full.table$generation = as.numeric(full.table$generation)
  full.table$str_generation = str_pad(full.table$generation, 5, pad = "0")
  full.table$name = paste(full.table$population, full.table$str_generation, full.table$strain, sep="-")
  full.table = full.table %>% mutate(strain = fct_reorder(strain, name))
  full.table$gene = factor(full.table$gene, levels=c("pykF", "nadR", "topA", "spoT", "ybaL", "hslU", "iclR", "yijC (fabR)"))

  ggplot(full.table, aes(x=gene, y=name, fill=mutated)) + geom_tile() + facet_grid(rows=vars(population), scales = "free")
  ggsave("LTEE_mutation_presence.pdf", width=8, height=24)

}
