#!/bin/bash
cd $1
#echo $1; \ # Debugging
module load samtools 2>/dev/null; \
find . -regex './*\.bam$' -exec sh -c " samtools view -H {} | grep '^@RG' |  tr '\n' '\t'; echo " \; |  sort | uniq | tr '\t' '\n'; \
find . -regex './*\.bam$' -exec sh -c " samtools flagstat {} | tr '\n' '\t'; echo " \; |  sort | uniq | tr '\t' '\n'