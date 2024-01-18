#!/bin/bash
set -e

source /opt/miniconda3/bin/activate base
mamba env update -f "${HOME}/repos/atac-seq/config/atac_env.yaml"
source /opt/miniconda3/bin/activate atac
source "${HOME}/repos/atac-seq/config/bash_src"
