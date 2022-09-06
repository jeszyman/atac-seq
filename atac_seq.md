
# Table of Contents

1.  [Prerequisites to run repository local integration testing](#orgd94d167)
2.  [Changelog](#org4a2beec)


<a id="orgd94d167"></a>

# Prerequisites to run repository local integration testing

-   Singularity container built from <https://github.com/jeszyman/atac-seq/blob/master/config/atac_Dockerfile>
-   Local snakemake
-   Local snakemake configuration YAML


<a id="org4a2beec"></a>

# Changelog

-   <span class="timestamp-wrapper"><span class="timestamp">[2022-09-06 Tue] </span></span> Re-written for my biotools repo best practices <span class="timestamp-wrapper"><span class="timestamp">[2022-09-06 Tue]</span></span>. Downgraded to alignment and qc only. Need to add back macs2.
-   <span class="timestamp-wrapper"><span class="timestamp">[2022-08-29 Mon] </span></span> Initial pre-processing, peak calling, and normalization validated.
