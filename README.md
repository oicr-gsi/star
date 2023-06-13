# star

STAR 2.1

## Overview

## Dependencies

* [star 2.7.10b](https://github.com/alexdobin/STAR)
* [picard 2.19.2](https://broadinstitute.github.io/picard/)


## Usage

### Cromwell
```
java -jar cromwell.jar run star.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputGroups`|Array[InputGroup]|Array of fastq files to align with STAR and the merged filename
`outputFileNamePrefix`|String|Prefix for filename
`reference`|String|Reference id, hg19 or hg38


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`runStar.genomeIndexDir`|String|"$HG38_STAR_INDEX100_ROOT/"|Path to STAR index
`runStar.starSuffix`|String|"Aligned.sortedByCoord.out"|Suffix for sorted file
`runStar.transcriptomeSuffix`|String|"Aligned.toTranscriptome.out"|Suffix for transcriptome-aligned file
`runStar.chimericjunctionSuffix`|String|"Chimeric.out"|Suffix for chimeric junction file
`runStar.genereadSuffix`|String|"ReadsPerGene.out"|ReadsPerGene file suffix
`runStar.addParam`|String?|None|Additional STAR parameters
`runStar.modules`|String|"hg38-star-index100/2.7.10b"|modules for running STAR
`runStar.chimOutType`|String|"WithinBAM SoftClip Junctions"|Indicate where chimeric reads are to be written
`runStar.outFilterMultimapNmax`|Int|50|max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
`runStar.chimScoreDropMax`|Int|30|max drop (difference) of chimeric score (the sum of scores of allchimeric segments) from the read length
`runStar.uniqMAPQ`|Int|255|Score for unique mappers
`runStar.saSparsed`|Int|2|saSparsed parameter for STAR
`runStar.multiMax`|Int|-1|multiMax parameter for STAR
`runStar.chimSegmin`|Int|10|minimum length of chimeric segment length
`runStar.chimJunOvMin`|Int|10|minimum overhang for a chimeric junction
`runStar.alignSJDBOvMin`|Int|10|minimum overhang for annotated spliced alignments
`runStar.alignMatGapMax`|Int|100000|maximum gap between two mates
`runStar.alignIntMax`|Int|100000|maximum intron size
`runStar.chimMulmapScoRan`|Int|3|the score range for multi-mapping chimeras below the best chimeric score
`runStar.chimScoJunNonGTAG`|Int|-1|penalty for a non-GTAG chimeric junction
`runStar.chimScoreSeparation`|Int|1|minimum difference (separation) between the best chimeric score and the next one
`runStar.chimMulmapNmax`|Int|50|maximum number of chimeric multi-alignments
`runStar.chimNonchimScoDMin`|Int|10|to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value
`runStar.chimOutJunForm`|Int?|None|flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata
`runStar.peOvNbasesMin`|Int|10|minimum number of overlap bases to trigger mates merging and realignment
`runStar.chimSegmentReadGapMax`|Int|3|maximum gap in the read sequence between chimeric segments
`runStar.peOvMMp`|Float|0.1|maximum proportion of mismatched bases in the overlap area
`runStar.threads`|Int|6|Requested CPU threads
`runStar.jobMemory`|Int|64|Memory allocated for this job
`runStar.timeout`|Int|72|hours before task timeout
`indexBam.jobMemory`|Int|12|Memory allocated indexing job
`indexBam.modules`|String|"picard/2.19.2"|modules for running indexing job
`indexBam.timeout`|Int|48|hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`starBam`|File|Output bam aligned to genome
`starIndex`|File|Output index file for bam aligned to genome
`transcriptomeBam`|File|Output bam aligned to transcriptome
`starChimeric`|File|Output chimeric junctions file
`geneReadFile`|File|Output raw read counts per transcript


## Commands
 This section lists command(s) run by the STAR workflow.
 
 ### run STAR aligner
 
  STAR --twopassMode Basic \
       --genomeDir ~{genomeIndexDir} \
       --readFilesIn ~{sep="," read1s} ~{sep="," read2s} \
       --readFilesCommand zcat \
       --outFilterIntronMotifs RemoveNoncanonical \
       --outFileNamePrefix ~{outputFileNamePrefix}. \
       --outSAMmultNmax ~{multiMax} \
       --outSAMattrRGline ~{sep=" , " readGroups} \
       --outSAMstrandField intronMotif \
       --outSAMmapqUnique  ~{uniqMAPQ} \
       --outSAMunmapped Within KeepPairs \
       --genomeSAsparseD ~{saSparsed} \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode TranscriptomeSAM GeneCounts \
       --chimSegmentMin ~{chimSegmin} \
       --chimJunctionOverhangMin ~{chimJunOvMin} \
       --alignSJDBoverhangMin ~{alignSJDBOvMin} \
       --alignMatesGapMax ~{alignMatGapMax} \
       --alignIntronMax ~{alignIntMax} \
       --alignSJstitchMismatchNmax 5 -1 5 5 \
       --chimMultimapScoreRange ~{chimMulmapScoRan} \
       --chimScoreJunctionNonGTAG ~{chimScoJunNonGTAG} \
       --chimMultimapNmax ~{chimMulmapNmax} \
       --chimNonchimScoreDropMin ~{chimNonchimScoDMin} \
       ~{"--chimOutJunctionFormat " + chimOutJunForm} \
       --peOverlapNbasesMin ~{peOvNbasesMin} \
       --peOverlapMMp ~{peOvMMp} \
       --outFilterMultimapNmax ~{outFilterMultimapNmax} \
       --runThreadN ~{threads} --chimOutType ~{chimOutType} \
       --chimScoreDropMax ~{chimScoreDropMax} \
       --chimScoreSeparation ~{chimScoreSeparation} \
       --chimSegmentReadGapMax ~{chimSegmentReadGapMax} ~{addParam}
 
 ### Process Chimeric junctions file for downstream use by STAR-Fusion
 
  awk 'NR<2{print $0;next}{print $0| "sort -V"}' ~{outputFileNamePrefix}.~{chimericjunctionSuffix}.junction \
  > tmp && mv tmp ~{outputFileNamePrefix}.~{chimericjunctionSuffix}.junction
 
 ### Index Bam file for random access
 
  java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar BuildBamIndex \
                               VALIDATION_STRINGENCY=LENIENT \
                               OUTPUT="~{basename(inputBam, '.bam')}.bai" \
                               INPUT=~{inputBam}
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
