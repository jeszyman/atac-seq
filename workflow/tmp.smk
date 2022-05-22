import pandas as pd            
import re

libraries = (
    pd.read_csv("/home/jeszyman/repos/atac-seq/test/inputs/full_libraries.tsv", sep="\t",
		dtype={"library_id": str})
    .set_index("library_id", drop=False)
    .sort_index()
)

FQ_DICT = dict(zip(libraries['library_id'], libraries['fq_basename']))

def fq_r1_path(library_id):
    return("/home/jeszyman/repos/atac-seq/test/inputs" + "/" + FQ_DICT[library_id] + "_R1.fastq.gz")

fq_path("lib001")

def fq_r1_path(library_id):
    return(config["fq_dir"] + "/" + FQ_DICT[library_id] + "_R1.fastq.gz")

container: config["container"]

FQ_BASENAME = pd.Series(libraries.fq_basename, dtype="string")

LIBRARY_IDS = ["atac1","atac2","atac3","atac4"]

MACS_BROAD_EXT = ["peaks.broadPeak", "peaks.gappedPeak", "peaks.xls"]

MACS_NARROW_EXT = ["peaks.narrowPeak", "summits.bed"]

BAM_PROCESS = ["regfilt", "open"]

rule all:
    input:
        expand(config["fq_dir"] + "/{library_id}_R1.fastq.gz", library_id = FQ_DICT.keys()),

def fq_r1_path(n):
    return("TEST" + lambda wcs: FQ_DICT[wcs])

test = fq_r1_path("lib001")

print(test("lib001"))

def myfunc(n):
  return lambda a : a * n

mydoubler = myfunc(2)

print(mydoubler(11))

rule symlink:
    params:
        r1 = lambda wcs: return("/home/jeszyman/repos/atac-seq/test/inputs" + "/" + FQ_DICT[wcs.f] + "_R1.fastq.gz"), 
    output:
        r1 = config["fq_dir"] + "/{f}_R1.fastq.gz",
    shell:
        """
        ln -sf {params.r1} {output.r1}
        """
