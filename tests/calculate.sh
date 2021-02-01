#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load samtools/1.9 2>/dev/null

ls | sort

find . -regex '.*\.bam$' -exec samtools view -H {} \; | grep '^@RG' | sort

find . -regex '.*\.bam$' -exec samtools flagstat {} \; | sort

find . -regex '.*\.tab$' -exec wc -l {} \;

find . -regex '.*\.junction$' -exec wc -l {} \;
