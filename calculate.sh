#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load samtools/1.9 2>/dev/null

ls | sort

find . -regex '.*\.bam$' -exec samtools view -H {} \; | grep '^@RG' | sort

find . -regex '.*\.bam$' -exec samtools flagstat {} \; | sort

find . -regex '.*\.bam$' -exec /bin/bash -c "samtools view {} | md5sum" \; | sort

find . -regex '.*\.bai$' -exec /bin/bash -c "samtools view {} | md5sum" \; | sort

find . -regex '.*\.Chimeric.out.junction$' -exec /bin/bash -c "cat {} | md5sum" \; | sort

find . -regex '.*\.ReadsPerGene.out.tab$' -exec /bin/bash -c "cat {} | md5sum" \; | sort
