#!/bin/bash
PROJECT_DIR="/home/ubuntu/bioinformatics_course"
RESULTS_DIR="${PROJECT_DIR}/results"
DATA_DIR="${PROJECT_DIR}/data"
ALIGN_DIR="${DATA_DIR}/aligned_data"
STATS_DIR="${DATA_DIR}/stats"

echo "emptying results"
cd $RESULTS_DIR
rm -r *
echo "emptying stats"
cd $STATS_DIR
rm -r *
echo "emptying aligned"
cd $ALIGN_DIR
rm -r *
echo "emptying trimmed"
cd $DATA_DIR/trimmed_fastq
rm -r *



