#!/usr/bin/env Rscript
#########1#########2#########3#########4#########5#########6#########7#########8
###
###    SCRIPT TITLE   ###
###

args = commandArgs(trailingOnly = TRUE)
passed_list = args[1]
saveloc = args[2]


filt_libs_raw = passed_list
filt_libs = unlist(strsplit(filt_libs_raw, " "))
saveloc = "/mnt/ris/jschwarz/cardiac-radiobiology/tmp/test.RData"

save(filt_libs_raw,
     file = saveloc)
