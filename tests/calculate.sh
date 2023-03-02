#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

module load samtools/1.9 2>/dev/null

ls | sort

find . -regex '.*\.bam$' -exec samtools view -H {} \; | grep '^@RG' | sort

find . -regex '.*\.bam$' -exec samtools flagstat {} \;
echo "md5sum for .tab file(s):"
find . -regex '.*\.tab$' -exec /bin/bash -c "cat {} | md5sum" \; | sort &&
echo "md5sum for .junction file(s):" &&
find . -regex '.*\.junction$' -exec /bin/bash -c "cat {} | grep -v '^#' | sort -V | md5sum" \; | sort
