args = commandArgs(trailingOnly = TRUE)
libraries_tsv = args[1]
groups_str = args[2]
macs2_dir = args[3]

library(tidyverse)

libraries = read_tsv(libraries_tsv)

groups = unlist(strsplit(groups_str, " "))

make_consensus = function(library_tib, group_filt){
  path_vect =
    library_tib %>%
    filter(group == group_filt) %>%
    mutate(path = paste0("test/macs2/", library, "_single_peaks.narrowPeak")) %>%
    pull(path)
  path_str = paste(path_vect, sep="", collapse=" ")
  system(
    paste0("bedops --intersect ", path_str, " > ",macs2_dir,"/",group_filt,"_consensus.bed")
    )
}

for (i in 1:length(groups)){
  make_consensus(libraries, groups[i])
}
