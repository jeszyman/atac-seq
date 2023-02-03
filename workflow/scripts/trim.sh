#!/usr/bin/env bash

# Functions
fastp_wrap(){
    local input_r1="${1}"
    local input_r2="${2}"
    local log_html="${3}"
    local log_json="${4}"
    local output_r1="${5}"
    local output_r2="${6}"
    local threads="${7}"
    #
    fastp --detect_adapter_for_pe \
          --html $log_html \
          --json $log_json \
          --in1 $input_r1 \
          --in2 $input_r2 \
          --length_limit 2000 \
          --out1 $output_r1 \
          --out2 $output_r2 \
          --thread $threads --trim_tail1 1
}

fastp_wrap $1 $2 $3 $4 $5 $6 $7
