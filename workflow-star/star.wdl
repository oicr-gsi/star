version 1.0

workflow starWorkflow {
input {
 File fastqR1
 File fastqR2
 String genomeIndexDir
 String? additionalParameters = ""
 Int? uniqMAPQ = 60
 Int? saSparsed = 2
 Int? multiMax = -1
 Int? threads = 4
 String RGID
 String RGLB
 String RGPL
 String RGPU
 String RGSM
 String? RGCM = "OICR"
}

call makeName{ input: fastqFile = fastqR1 }
call makeRg{ input: RGID = RGID, RGLB = RGLB, RGPL = RGPL, RGPU = RGPU, RGSM = RGSM, RGCM = RGCM }
call runStar{ input: fastqR1 = fastqR1, fastqR2 = fastqR2, genome_index_dir = genomeIndexDir, fileNAME = makeName.outputName, addParam  = additionalParameters, uniqMAPQ = uniqMAPQ, saSparsed = saSparsed, multiMax = multiMax, threads = threads, rgLine = makeRg.rgLine }
call indexBam as finalIndex{ input: inputBam = runStar.outputBam }


output {
  File starBam   = runStar.outputBam
  File starIndex = finalIndex.outputBai
  File transcriptomeBam = runStar.transcriptomeBam
  File geneReadFile     = runStar.geneReads
}
}

# ==========================================
#   MAKE NAME
# ==========================================
task makeName {
input {
 File fastqFile
}

command <<<
 FILE_NAME=$(echo ~{basename(fastqFile, ".fastq.gz")} | sed s/_R1.*// )
 echo $FILE_NAME 
>>>

output {
 String outputName = read_string(stdout())
}
}

# ==========================================
#   MAKE RGLINE
# ==========================================
task makeRg {
input {
 String RGID
 String RGLB
 String RGPL
 String RGPU
 String RGSM
 String? RGCM
}

command <<<
 RG=$(echo "ID:~{RGID} PL:~{RGPL} PU:~{RGPU} LB:~{RGLB} SM:~{RGSM} CM:~{RGCM}")
 echo $RG
>>>

output {
 String rgLine = read_string(stdout())
}
}


# ==========================================
#  TASK 1 of 2: run STAR aligner
# ==========================================
task runStar {
input {
  File fastqR1
  File fastqR2
  String genome_index_dir
  String fileNAME
  String rgLine
  String? starSuffix = "Aligned.sortedByCoord.out"
  String? transcriptomeSuffix = "Aligned.toTranscriptome.out"
  String? genereadSuffix = "ReadsPerGene.out"
  String? addParam = ""
  String? modules = "star/2.6.0c"
  Int? uniqMAPQ = 60
  Int? saSparsed = 2
  Int? multiMax = -1
  Int? threads = 6
  Int? jobMemory  = 36
}
# missing --clip3pAdapterSeq $adaptors
command <<<
 $STAR_ROOT/bin/STAR \
      --twopassMode Basic \
      --genomeDir ~{genome_index_dir} \
      --readFilesIn ~{fastqR1} ~{fastqR2} \
      --readFilesCommand zcat \
      --outFilterIntronMotifs RemoveNoncanonical \
      --outFileNamePrefix ~{fileNAME} \
      --outSAMmultNmax ~{multiMax} \
      --outSAMattrRGline ~{rgLine} \
      --outSAMmapqUnique  ~{uniqMAPQ} \
      --outSAMunmapped Within KeepPairs \
      --genomeSAsparseD "~{saSparsed}" \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM GeneCounts \
      --runThreadN "~{threads}" ~{addParam}
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  cpu:     "~{threads}"
}


output {
 File outputBam        = "~{fileNAME}~{starSuffix}.bam"
 File transcriptomeBam = "~{fileNAME}~{transcriptomeSuffix}.bam"
 File geneReads        = "~{fileNAME}~{genereadSuffix}.tab"
}
}

# ==========================================
#  TASK 2 of 2: index bam file with picard
# ==========================================
task indexBam {
input {
	File   inputBam
        Int?   jobMemory  = 8
        Int?   javaMemory = 6
        String? modules = "java/8 picard/2.19.2" 
}

command <<<
 java -Xmx~{javaMemory}G -jar $PICARD_ROOT/picard.jar BuildBamIndex \
                              VALIDATION_STRINGENCY=LENIENT \
                              OUTPUT="~{basename(inputBam, '.bam')}.bai" \
                              INPUT=~{inputBam} 
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
  File outputBai = "~{basename(inputBam, '.bam')}.bai"
}
}

