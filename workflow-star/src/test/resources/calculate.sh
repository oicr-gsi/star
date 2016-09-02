#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

qrsh -l h_vmem=4G -cwd -now n "\
. /oicr/local/Modules/default/init/bash ; \
module load samtools 2>/dev/null; \
find . -regex './data/.*.bam$' -exec sh -c \" samtools view -H {} | grep '^@RG' |  tr '\n' '\t'; echo \" \; \
| sort | uniq | tr '\t' '\n'"
qrsh -l h_vmem=4G -cwd -now n "\
. /oicr/local/Modules/default/init/bash ; \
module load samtools 2>/dev/null; \
find . -regex './data/.*.bam$' -exec sh -c \" samtools flagstat {} | tr '\n' '\t'; echo \" \; \
| sort | uniq | tr '\t' '\n'"