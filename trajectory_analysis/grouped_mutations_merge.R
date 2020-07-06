grouped_mutations_merge <- function(input_path, input_file_name, grouped_mutations_file_name, output_path) {

source("common.R")
dir.create(file.path(output_path), showWarnings = FALSE)
  
cat(input_path, input_file_name, grouped_mutations_file_name, output_path, "\n")

read_counts = read.table(file.path(input_path, input_file_name), header=T)
grouped_mutations = read.csv(grouped_mutations_file_name, header=T)
grouped_mutations$full_name = as.character(grouped_mutations$full_name)
grouped_mutations$resolved = as.character(grouped_mutations$resolved)

ncol(read_counts)

for (i in 1:nrow(grouped_mutations)) {
  matched_col_index = match(grouped_mutations$full_name[i], colnames(read_counts))
  if (!is.na(matched_col_index)) {
    cat("Handling... ", grouped_mutations$full_name[i], "\n")
    if (grouped_mutations$resolved[i] == "merged") {
      cat("Deleting b/c merged: ", colnames(read_counts)[existing_col_index], "\n")
      read_counts = read_counts[-matched_col_index]
    } else {
      existing_col_index = match(grouped_mutations$resolved[i], colnames(read_counts))
      if (!is.na(existing_col_index) && (existing_col_index != matched_col_index) ) {
        cat("Merging column:",  colnames(read_counts)[matched_col_index], " to ", colnames(read_counts)[existing_col_index], "\n")
        # the merged name already exists... add together the column numbers and delete this one
        cat(read_counts[,existing_col_index], "\n")
        read_counts[,existing_col_index] = read_counts[,existing_col_index] + read_counts[,matched_col_index]
        cat(read_counts[,existing_col_index], "\n")
        read_counts = read_counts[-matched_col_index]
      } else {
        cat("Renamed column: ", colnames(read_counts)[matched_col_index], " to ", grouped_mutations$resolved[i], "\n")
        colnames(read_counts)[matched_col_index] = grouped_mutations$resolved[i]
      }
    }
  }
}

ncol(read_counts)
write.table(read_counts, file.path(output_path, input_file_name), sep="\t", row.names=F)
}
