
## - Rscript

args = commandArgs(trailingOnly = TRUE)
bam_list_str = args[1]
peak_ratio_tsv = args[2]
peak_ratio_plot = args[3]
peak_cut = args[4]

library(GenomicAlignments)
library(tidyverse)

bam_list = unlist(strsplit(bam_list_str, " "))
names(bam_list) = substr(gsub("^.*lib","lib", bam_list), 1, 6)

tally_lengths = function(in_bam){
  # Make a tibble with counts of fragment lengths
  gal = readGAlignments(in_bam,
                        param=ScanBamParam(what=c("isize")))
  tib = mcols(gal) %>%
    as_tibble() %>%
    mutate(frag_len = abs(isize)) %>%
    group_by(frag_len) %>%
    tally()
  return(tib)
}

change_column_name <- function(x, aList) {
  # For each per-library tibble, change column to library ID
  dat <- aList[[x]]
  names(dat)[2] <- x
  return(dat)
}

pre_frag_len_list = lapply(bam_list, tally_lengths)
frag_len_list = lapply(names(pre_frag_len_list), change_column_name, pre_frag_len_list)

frags =
  # Make a complete list of fragment sizes, 1-1000
  data.frame(frag_len = 1:1000) %>%
  as_tibble()
frags

frag_len_tib =
  frag_len_list %>% purrr::reduce(full_join, by = "frag_len") %>%
  full_join(frags, by = "frag_len") %>%
  arrange(frag_len) %>%
   replace(is.na(.), 0) %>%
   mutate(frag_len_fct = ifelse(frag_len > 1000, "other", frag_len)) %>%
  select(!frag_len) %>%
  pivot_longer(!frag_len_fct, names_to = "library", values_to = "count") %>%
  group_by(frag_len_fct, library) %>%
  summarize(count = sum(count)) %>%
  mutate(frag_len_fct = as.numeric(frag_len_fct)) %>%
  arrange(frag_len_fct)

cut = frag_len_tib %>%
  mutate(mono_cut = ifelse(frag_len_fct < 146, "open", "mono")) %>%
  group_by(library,mono_cut) %>%
  summarize(high = max(count)) %>%
  pivot_wider(names_from = mono_cut, values_from = high) %>%
  mutate(peak_ratio = open / mono) %>%
  select(library, open, mono, peak_ratio)

cut %>% write_tsv(file = peak_ratio_tsv)

plot = frag_len_tib %>%
  left_join(cut, by = "library") %>%
  mutate(peak_ratio_mod = ifelse(peak_ratio < peak_cut, NA, peak_ratio)) %>%
  ggplot(., aes(x=frag_len_fct, y = count)) +
  geom_line(aes(color = peak_ratio_mod)) +
  facet_wrap(vars(library)) +
  xlab("Fragment Length") + ylab("Count") +
  geom_hline(aes(yintercept = mono)) +
  scale_color_continuous(name = "Ratio of open to mononucleosomal peaks", na.value = "red") +
  theme(legend.position = "bottom")

ggsave(plot, file = peak_ratio_plot)


## - Rscript

args = commandArgs(trailingOnly = TRUE)
bam_list_str = args[1]
peak_ratio_tsv = args[2]
peak_ratio_plot = args[3]
peak_cut = args[4]

library(GenomicAlignments)
library(tidyverse)

bam_list = unlist(strsplit(bam_list_str, " "))
names(bam_list) = substr(gsub("^.*lib","lib", bam_list), 1, 6)

tally_lengths = function(in_bam){
  # Make a tibble with counts of fragment lengths
  gal = readGAlignments(in_bam,
                        param=ScanBamParam(what=c("isize")))
  tib = mcols(gal) %>%
    as_tibble() %>%
    mutate(frag_len = abs(isize)) %>%
    group_by(frag_len) %>%
    tally()
  return(tib)
}

change_column_name <- function(x, aList) {
  # For each per-library tibble, change column to library ID
  dat <- aList[[x]]
  names(dat)[2] <- x
  return(dat)
}

pre_frag_len_list = lapply(bam_list, tally_lengths)
frag_len_list = lapply(names(pre_frag_len_list), change_column_name, pre_frag_len_list)

frags =
  # Make a complete list of fragment sizes, 1-1000
  data.frame(frag_len = 1:1000) %>%
  as_tibble()
frags

frag_len_tib =
  frag_len_list %>% purrr::reduce(full_join, by = "frag_len") %>%
  full_join(frags, by = "frag_len") %>%
  arrange(frag_len) %>%
   replace(is.na(.), 0) %>%
   mutate(frag_len_fct = ifelse(frag_len > 1000, "other", frag_len)) %>%
  select(!frag_len) %>%
  pivot_longer(!frag_len_fct, names_to = "library", values_to = "count") %>%
  group_by(frag_len_fct, library) %>%
  summarize(count = sum(count)) %>%
  mutate(frag_len_fct = as.numeric(frag_len_fct)) %>%
  arrange(frag_len_fct)

cut = frag_len_tib %>%
  mutate(mono_cut = ifelse(frag_len_fct < 146, "open", "mono")) %>%
  group_by(library,mono_cut) %>%
  summarize(high = max(count)) %>%
  pivot_wider(names_from = mono_cut, values_from = high) %>%
  mutate(peak_ratio = open / mono) %>%
  select(library, open, mono, peak_ratio)

cut %>% write_tsv(file = peak_ratio_tsv)

plot = frag_len_tib %>%
  left_join(cut, by = "library") %>%
  mutate(peak_ratio_mod = ifelse(peak_ratio < peak_cut, NA, peak_ratio)) %>%
  ggplot(., aes(x=frag_len_fct, y = count)) +
  geom_line(aes(color = peak_ratio_mod)) +
  facet_wrap(vars(library)) +
  xlab("Fragment Length") + ylab("Count") +
  geom_hline(aes(yintercept = mono)) +
  scale_color_continuous(name = "Ratio of open to mononucleosomal peaks", na.value = "red") +
  theme(legend.position = "bottom")

ggsave(plot, file = peak_ratio_plot)
