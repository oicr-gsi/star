#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1
module load samtools/0.1.19 2>/dev/null
qrsh -V -l h_vmem=8G -cwd -now n "\
find . -regex '.*\.bam$' \
       -exec sh -c \" samtools flagstat {} | tr '\n' '\t'; echo \" \; \
| sort | uniq | tr '\t' '\n'"
