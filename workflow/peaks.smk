print("Integration testing snakefile for ATAC-seq peak analysis\n")

# Import common packages
import pandas as pd
import re
import numpy as np





rule all:
    input:

rule symlink_peaks_inputs:
    input:

include: config["atac_repo"] + "/workflow/atac.smk"
