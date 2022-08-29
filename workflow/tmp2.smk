import pandas as pd            

libraries = (
    pd.read_csv("/home/jeszyman/repos/atac-seq/test/inputs/full_libraries.tsv", sep="\t",
		dtype={"library_id": str})
    .set_index("library_id", drop=False)
    .sort_index()
)

rule all:
    input:
       expand(config["fq_sym_dir"] + "/{library_id}_{read}.fastq.gz", library_id = libraries.library_id, read = ["R1", "R2"])

       
rule symlink_fastqs:
    params:
        fastq = lambda w: libraries[libraries.library_id == w.library_id].fq_basename.tolist()
    output:
        r1 = config["fq_sym_dir"] + "/{library_id}_R1.fastq.gz",
        r2 = config["fq_sym_dir"] + "/{library_id}_R2.fastq.gz",	
    shell:
        """
        ln -sf --relative {config[fq_src_dir]}/{params.fastq}_R1.fastq.gz {output.r1}
        ln -sf --relative {config[fq_src_dir]}/{params.fastq}_R2.fastq.gz {output.r2}
        """

