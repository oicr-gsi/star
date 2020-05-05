version 1.0

workflow star {
input {
 Array[Pair[Pair[File, File], String]]+ inputFqsRgs
 String outputFileNamePrefix
}

scatter (FqRg in inputFqsRgs) {
    File read1s       = FqRg.left.left
    File read2s       = FqRg.left.right
    String readGroups = FqRg.right
}

call runStar { input: read1s = read1s, read2s = read2s, readGroups = readGroups, outputFileNamePrefix = outputFileNamePrefix }
call indexBam { input: inputBam = runStar.outputBam }

meta {
 author: "Peter Ruzanov"
 email: "peter.ruzanov@oicr.on.ca"
 description: "STAR 2.0"
 dependencies: [
      {
        name: "star/2.7.3a",
        url: "https://github.com/alexdobin/STAR"
      },
      {
        name: "picard/2.19.2",
        url: "https://broadinstitute.github.io/picard/"
      }
    ]
}

output {
  File starBam          = runStar.outputBam
  File starChimeric     = runStar.outputChimeric
  File starIndex        = indexBam.outputBai
  File transcriptomeBam = runStar.transcriptomeBam
  File geneReadFile     = runStar.geneReads
 }
}

# ==========================================
#  TASK 1 of 2: run STAR aligner
# ==========================================
task runStar {
input {
  Array[File]+ read1s
  Array[File]+ read2s
  Array[String]+ readGroups
  String genomeIndexDir = "$HG38_STAR_INDEX100_ROOT/"
  String outputFileNamePrefix
  String starSuffix = "Aligned.sortedByCoord.out"
  String transcriptomeSuffix = "Aligned.toTranscriptome.out"
  String chimericjunctionSuffix = "Chimeric.out"
  String genereadSuffix = "ReadsPerGene.out"
  String? addParam
  String modules = "star/2.7.3a hg38-star-index100/2.7.3a"
  Int uniqMAPQ = 255
  Int saSparsed = 2
  Int multiMax = -1
  Int chimSegmin = 12
  Int chimJunOvMin = 12
  Int alignSJDBOvMin = 10
  Int alignMatGapMax = 100000
  Int alignIntMax = 100000
  Int chimMulmapScoRan = 3
  Int chimScoJunNonGTAG = -4
  Int chimMulmapNmax = 20
  Int chimNonchimScoDMin = 10
  Int? chimOutJunForm
  Int peOvNbasesMin = 12
  Float peOvMMp = 0.1
  Int threads = 6
  Int jobMemory = 64
  Int timeout = 72
}

parameter_meta {
 read1s: "array of read1s"
 read2s: "array of read2s"
 readGroups: "array of readgroup lines"
 starSuffix: "Suffix for sorted file"
 transcriptomeSuffix: "Suffix for transcriptome-aligned file"
 chimericjunctionSuffix: "Suffix for chimeric junction file"
 genereadSuffix: "ReadsPerGene file suffix"
 addParam: "Additional STAR parameters"
 modules: "modules for running STAR"
 uniqMAPQ: "Score for unique mappers"
 saSparsed: "saSparsed parameter for STAR"
 multiMax: "multiMax parameter for STAR"
 chimSegmin: "minimum length of chimeric segment length"
 chimJunOvMin: "minimum overhang for a chimeric junction"
 alignSJDBOvMin: "minimum overhang for annotated spliced alignments"
 alignMatGapMax: "maximum gap between two mates"
 alignIntMax: "maximum intron size"
 chimMulmapScoRan: "the score range for multi-mapping chimeras below the best chimeric score"
 chimScoJunNonGTAG: "penalty for a non-GTAG chimeric junction"
 chimMulmapNmax: "maximum number of chimeric multi-alignments"
 chimNonchimScoDMin: "to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value"
 chimOutJunForm: "flag to add metadata to chimeric junction output for functionality with starFusion - 1 for metadata, 0 for no metadata"
 peOvNbasesMin: "minimum number of overlap bases to trigger mates merging and realignment"
 peOvMMp: "maximum proportion of mismatched bases in the overlap area"
 threads: "Requested CPU threads"
 jobMemory: "Memory allocated for this job"
 timeout: "hours before task timeout"
}

# missing --clip3pAdapterSeq $adaptors
command <<<
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
      --runThreadN ~{threads} ~{addParam}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
  timeout: "~{timeout}"
}

output {
 File outputBam        = "~{outputFileNamePrefix}.~{starSuffix}.bam"
 File outputChimeric   = "~{outputFileNamePrefix}.~{chimericjunctionSuffix}.junction"
 File transcriptomeBam = "~{outputFileNamePrefix}.~{transcriptomeSuffix}.bam"
 File geneReads        = "~{outputFileNamePrefix}.~{genereadSuffix}.tab"
}

meta {
  output_meta: {
    outputBam:        "Output bam aligned to genome",
    outputChimeric:   "Output chimeric junctions file",
    transcriptomeBam: "Output bam aligned to transcriptome",
    geneReads:        "Output raw read counts per transcript"
  }
}

}

# ==========================================
#  TASK 2 of 2: index bam file with picard
# ==========================================
task indexBam {
input {
	File  inputBam
  Int   jobMemory = 12
  String modules  = "picard/2.19.2"
  Int timeout     = 48
}

parameter_meta {
 inputBam:  "Input bam file"
 jobMemory: "Memory allocated indexing job"
 modules:   "modules for running indexing job"
 timeout:   "hours before task timeout"
}

command <<<
 java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar BuildBamIndex \
                              VALIDATION_STRINGENCY=LENIENT \
                              OUTPUT="~{basename(inputBam, '.bam')}.bai" \
                              INPUT=~{inputBam} 
>>>

runtime {
   memory: "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File outputBai = "~{basename(inputBam, '.bam')}.bai"
}

meta {
  output_meta: {
    outputBai: "Output index file for bam aligned to genome"
  }
}

}