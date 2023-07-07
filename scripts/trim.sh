#!/usr/bin/env bash

input_r1="${1}"
input_r2="${2}"
log_html="${3}"
log_json="${4}"
output_r1="${5}"
output_r2="${6}"
threads="${7}"

# Functions
fastp_wrap(){
    #
    fastp --detect_adapter_for_pe \
          --html $log_html \
          --json $log_json \
          --in1 $input_r1 \
          --in2 $input_r2 \
          --out1 $output_r1 \
          --out2 $output_r2 \
          --thread $threads --trim_tail1 1
}

fastp_wrap $input_r1 \
           $input_r2 \
           $log_html \
           $log_json \
           $output_r1 \
           $output_r2 \
           $threads
