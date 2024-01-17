#!/bin/bash
source /opt/miniconda3/bin/activate base
mamba env update config/atac_env.yaml
source /opt/miniconda3/bin/activate atac
source ./config/bash_src
