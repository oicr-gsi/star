version 1.0

workflow star {
input {
 File fastqR1
 File? fastqR2
 String? additionalParameters = ""
 String? outputFileNamePrefix = ""
 String RGID
 String RGLB
 String RGPL
 String RGPU
 String RGSM
 String? RGCM
}

String? outputPrefix = if outputFileNamePrefix=="" then basename(fastqR1, '.fastq.gz') else outputFileNamePrefix

call makeRg{ input: RGID = RGID, RGLB = RGLB, RGPL = RGPL, RGPU = RGPU, RGSM = RGSM, RGCM = RGCM }
call runStar{ input: fastqR1 = fastqR1, fastqR2 = fastqR2, fileNAME = outputPrefix, addParam  = additionalParameters, rgLine = makeRg.rgLine }
call indexBam as finalIndex{ input: inputBam = runStar.outputBam }

meta {
 author: "Peter Ruzanov"
 email: "peter.ruzanov@oicr.on.ca"
 description: "STAR 2.0"
}

output {
  File starBam   = runStar.outputBam
  File starIndex = finalIndex.outputBai
  File transcriptomeBam = runStar.transcriptomeBam
  File geneReadFile     = runStar.geneReads
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
 String? RGCM="MyCompany"
}

parameter_meta {
 RGID: "ID field in Read Group field"
 RGLB: "LB field in Read Group field"
 RGPL: "PL field in Read Group field"
 RGPU: "PU field in Read Group field"
 RGSM: "SM field in Read Group field"
 RGCM: "CM field in Read Group field"
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
  File? fastqR2
  String? genome_index_dir = "$HG19_STAR_INDEX100_ROOT/"
  String? fileNAME
  String rgLine
  String? starSuffix = "Aligned.sortedByCoord.out"
  String? transcriptomeSuffix = "Aligned.toTranscriptome.out"
  String? genereadSuffix = "ReadsPerGene.out"
  String? addParam = ""
  String? modules = "star/2.6.0c hg19-star-index100/2.6.0c"
  Int? uniqMAPQ = 60
  Int? saSparsed = 2
  Int? multiMax = -1
  Int? threads = 6
  Int? jobMemory  = 36
}

parameter_meta {
 fastqR1: "Mate 1 reads in fastq file"
 fastqR2: "Mate 2 reads in fastq file"
 genome_index_dir: "Directory with indices for STAR"
 fileNAME: "Prefix for building output file name"
 rgLine: "Line for making ReadGroup field"
 starSuffix: "Suffix for sorted file"
 transcriptomeSuffix: "Suffix for transcriptome-aligned file"
 genereadSuffix: "ReadsPerGene file suffix"
 addParam: "Additional STAR parameters"
 modules: "modules for running STAR"
 uniqMAPQ: "Score for unique mappers"
 saSparsed: "saSparsed parameter for STAR"
 multiMax: "multiMax parameter for STAR"
 threads: "Requested CPU threads"
 jobMemory: "Memory allocated for this job"
}

# missing --clip3pAdapterSeq $adaptors
command <<<
 STAR --twopassMode Basic \
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
        Int?   jobMemory  = 12
        String? modules = "java/8 picard/2.19.2" 
}

parameter_meta {
 inputBam: "Input bam file"
 jobMemory: "Memory allocated indexing job"
 modules: "modules for running indexing job"
}

command <<<
 java -Xmx~{jobMemory-6}G -jar $PICARD_ROOT/picard.jar BuildBamIndex \
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

