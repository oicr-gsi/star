package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.SemanticWorkflow;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class STARWorkflow extends SemanticWorkflow {

    String java;
    String picard_dir;
    String input1_path = null;
    String input2_path = null;
    String genome_index_dir = null;
    String dataDir = null;
    String outputFileName = null;
    boolean manualOutput;
    boolean produceTranscriptomeBam;
    String read1_adapterTrim = null;
    String read2_adapterTrim = null;
    String trimmedFile_1;
    String trimmedFile_2;
    //RG and STAR parameters
    String star;
    String RGID;
    String RGLB;
    String RGPL;
    ArrayList<String> RGPU = new ArrayList<String>();
    String RGSM;
    String RGCM;
    String additionalStarParams;
    int readTrimming; //aln
    int numOfThreads; //aln 
    String star_aln_params;
    ArrayList<String> read1 = new ArrayList<String>();
    ArrayList<String> read2 = new ArrayList<String>();
//    Boolean runMultipleFiles = false;
    SqwFile outputFile;
    SqwFile outputFileIndex;

    String queue;
    private int uniqMAPQ;
    private int multiMax;
    private int saSparsed;
    private final static String DEFAULT_UMAPQ = "60"; // For compatibility with downstream tools such as GATK
    private final static String DEFAULT_MULTI = "-1";
    private final static String DEFAULT_SASPARSED = "2";
    private final static String DEFAULT_THREADS = "6";
    private final static String STAR_SUFFIX = "Aligned.sortedByCoord.out";
    private final static String TRANSCRIPTOME_SUFFIX = "Aligned.toTranscriptome.out";
    private final static String GENEREAD_SUFFIX = "ReadsPerGene.out";
    //Ontology-related variables
    private static final String EDAM = "EDAM";
    private static final Map<String, Set<String>> cvTerms;
    
        static {
        cvTerms = new HashMap<String, Set<String>>();
        cvTerms.put(EDAM, new HashSet<String>(Arrays.asList("BAM", "BAI","Text","Alignment format",
                                                            "Sequence alignment","Gene expression")));
    }
        
    /**
     * Function that returns CV terms put into a Map container
     * @return myTerms
     */
    @Override
    protected Map<String, Set<String>> getTerms() {
       Map<String, Set<String>> myTerms = new HashMap<String, Set<String>>();
       myTerms.putAll(cvTerms);
       return myTerms;
    }
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {

            input1_path = getProperty("input_file_1");
            input2_path = getProperty("input_file_2");
            genome_index_dir = getProperty("index_dir");
            java = getProperty("bundled_jre");
            picard_dir = getProperty("picard_dir");
            if (!picard_dir.endsWith("/")) {picard_dir+="/";}
            star = getProperty("star");
           
            manualOutput = Boolean.valueOf(getOptionalProperty("manual_output", "false"));
            produceTranscriptomeBam = Boolean.valueOf(getOptionalProperty("produce_transcriptome_bam","true"));
            numOfThreads = Integer.valueOf(getOptionalProperty("star_aln_threads", DEFAULT_THREADS));
            uniqMAPQ     = Integer.valueOf(getOptionalProperty("uniqMapQ", DEFAULT_UMAPQ));
            multiMax     = Integer.valueOf(getOptionalProperty("multimap_max", DEFAULT_MULTI));
            saSparsed    = Integer.valueOf(getOptionalProperty("sa_sparsed", DEFAULT_SASPARSED));

            
            RGID = getProperty("rg_platform_unit");
            if (RGID.contains(" ")) {RGID = "\"" + RGID + "\"";}
            RGLB = getProperty("rg_library");
            if (RGLB.contains(" ")) {RGLB = "\"" + RGLB + "\"";}
            RGPL = getProperty("rg_platform");
            if (RGPL.contains(" ")) {RGPL = "\"" + RGPL + "\"";}
//            RGPU = getProperty("rg_platform_unit");
            String runDetails = getProperty("rg_platform_unit");
            if (runDetails.contains(" ")) {runDetails = "\"" + runDetails + "\"";}
            if (runDetails.contains(",")) {
                for (String rgpu:runDetails.split(",")){
                    RGPU.add(rgpu);
                }
            }
            else {
                RGPU.add(runDetails);
            }
//            if (RGPU.contains(" ")) {RGPU = "\"" + RGPU + "\"";}
            // split RGPU by sample entries
            RGSM = getProperty("rg_sample_name");
            if (RGSM.contains(" ")) {RGSM = "\"" + RGSM + "\"";}
            RGCM = getProperty("rg_organization");
            if (RGCM.contains(" ")) {RGCM = "\"" + RGCM + "\"";}
            additionalStarParams   = getOptionalProperty("additionalStarParams", "");

            queue = getOptionalProperty("queue", "");

            if (hasPropertyAndNotNull("outputFileName") && !getProperty("outputFileName").isEmpty()) {
                outputFileName = getProperty("outputFileName");
            } else {
		outputFileName = "SWID_" + getProperty("ius_accession") + "_" 
             	+ getProperty("rg_library") + "_" + getProperty("sequencer_run_name") + "_" + getProperty("barcode") 
             	+ "_L00" + getProperty("lane");
            }


        } catch (Exception e) {
            Logger.getLogger(STARWorkflow.class.getName()).log(Level.SEVERE, null, e);
            System.exit(1);
        }

        //TODO in a future we may accept multiple fastq and RG params (multi-lane alignment)
        // registers the first input file
       
        if (input1_path.contains(",")){
            String[] pathDesc = input1_path.split(",");
//            String[] provisionedPaths1 = new String[pathDesc.length];
            for (int id=0; id < pathDesc.length; id++ ){  
                String inPth = pathDesc[id];
                SqwFile r1 = this.createFile("file_in_R1" + id);
                r1.setSourcePath(inPth);
                r1.setType("chemical/seq-na-fastq-gzip");
                r1.setIsInput(true);
//                provisionedPaths1[id] = r1.getProvisionedPath();
                read1.add(r1.getProvisionedPath());
            }
        }
        else {
            SqwFile r1 = this.createFile("file_in_0");
            r1.setSourcePath(input1_path);
            r1.setType("chemical/seq-na-fastq-gzip");
            r1.setIsInput(true);
//            String[] provisionedPaths1 = new String[]{r1.getProvisionedPath()};
            read1.add(r1.getProvisionedPath());
        }
//        this.read1 = pr
        
        if (input2_path.contains(",")){
            String[] pathDesc = input2_path.split(",");
//            String[] provisionedPaths2 = new String[pathDesc.length];
            for (int id=0; id < pathDesc.length; id++ ){  
                String inPth = pathDesc[id];
                SqwFile r2 = this.createFile("file_in_R2" + id);
                r2.setSourcePath(inPth);
                r2.setType("chemical/seq-na-fastq-gzip");
                r2.setIsInput(true);
//                provisionedPaths2[id] = r2.getProvisionedPath();
                read2.add(r2.getProvisionedPath());
            }
        }
        else{
            // registers the second input file
            SqwFile r2= this.createFile("file_in_1");
            r2.setSourcePath(input2_path);
            r2.setType("chemical/seq-na-fastq-gzip");
            r2.setIsInput(true);
//            String[] provisionedPaths2 = new String[]{r2.getProvisionedPath()};
            read2.add(r2.getProvisionedPath());
        }
        
        

        outputFile = createOutputFile(this.dataDir + outputFileName + "." + STAR_SUFFIX + ".bam", "application/bam", manualOutput);
	outputFileIndex = createOutputFile(this.dataDir + outputFileName + "." + STAR_SUFFIX + ".bai", "application/bam-index", manualOutput);

        return this.getFiles();
    }

    @Override
    public void setupDirectory() {
        dataDir = getOptionalProperty("data_dir","data");
        if (!dataDir.endsWith("/")) {dataDir+="/";}
        this.addDirectory(dataDir);
    }

    @Override
    public void buildWorkflow() {
        
//        StringBuilder sb1 = new StringBuilder();
//        for (String r:read1){
//            sb1.append(r);
//        }
//        StringBuilder sb2 = new StringBuilder();
//        for (String r:read2){
//            sb2.append(r);
//        }
        String r1 = String.join(",", read1);
        String r2 = String.join(",", read2);
//        String r1 = sb1.toString();
//        String r2 = sb2.toString();       
        Job job01 = this.getWorkflow().createBashJob("star_align");

        job01.setCommand(star
                + " --twopassMode Basic"
                + " --genomeDir "   + this.genome_index_dir
                + " --readFilesIn " + r1 + " " + r2
                + " --readFilesCommand zcat "
                + " --outFilterIntronMotifs RemoveNoncanonical "
                + " --outFileNamePrefix " + this.dataDir + outputFileName + "."
                + " --outSAMmultNmax " + this.multiMax
                + " --outSAMattrRGline" + this.prepareRGLine() // create multiple --outSAMattrRGline ID:SAMPLE1 CN:XX DS:XXX , ID:SAMPLE2 CN:XX DS:XXX
                + " --outSAMmapqUnique " + this.uniqMAPQ
                + " --outSAMunmapped Within KeepPairs "
                + " --genomeSAsparseD " + this.saSparsed
                + " --outSAMtype BAM SortedByCoordinate");
        String addParams = parameters();
        if (null != addParams) {
            job01.getCommand().addArgument(addParams);
        }
        job01.setMaxMemory(getProperty("star_aln_mem_mb"));
        job01.setQueue(queue);
        
        if (this.numOfThreads > 1) {
         job01.getCommand().addArgument(" --runThreadN " + this.numOfThreads);
         job01.setThreads(this.numOfThreads);
        }
        
        this.attachCVterms(outputFile, EDAM, "BAM,Sequence alignment,Alignment format");
	job01.addFile(outputFile);
        
        Job job02 = this.indexBamJob(this.dataDir + outputFileName + "." + STAR_SUFFIX);
        job02.addParent(job01);
        
        this.attachCVterms(outputFileIndex, EDAM, "BAI,Sequence alignment,Alignment format");
	job02.addFile(outputFileIndex);
        
        if (this.produceTranscriptomeBam) {
            
            SqwFile transBamFile = createOutputFile(this.dataDir + outputFileName + "." + TRANSCRIPTOME_SUFFIX + ".bam", "application/bam", manualOutput);
            this.attachCVterms(transBamFile, EDAM, "BAM,Sequence alignment,Alignment format");
	    job01.addFile(transBamFile);
            SqwFile geneReadFile = createOutputFile(this.dataDir + outputFileName + "." + GENEREAD_SUFFIX + ".tab", "text/plain", manualOutput);
            this.attachCVterms(geneReadFile, EDAM, "Text,Gene expression");
	    job01.addFile(geneReadFile);
            
        }

    }
    
    public String prepareRGLine() {
        StringBuilder sb = new StringBuilder();
        for (String rgpu:RGPU){
//            String[] rgBuild = new String[]() {"ID:" + RGID};
            sb.append("ID:").append(RGID).append(" ");
            sb.append("PL:").append(RGPL).append(" ");
            sb.append("PU:").append(rgpu).append(" ");
            sb.append("LB:").append(RGLB).append(" ");
            sb.append("SM:").append(RGSM).append(" ");
            sb.append("CM:").append(RGCM).append(" ");
            sb.append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
    
    /**
     * Note that file path should not have extension
     * @param inputFile
     * @return 
     */
    protected Job indexBamJob(String inputFile) {
    Job jobIndex = this.getWorkflow().createBashJob("index_bam");
        jobIndex.setCommand(this.java + " -Xmx3G -jar "
                + this.picard_dir + "BuildBamIndex.jar"
                + " I=" + inputFile + ".bam"
                + " O=" + inputFile + ".bai"
                + " VALIDATION_STRINGENCY=LENIENT"); // The last one is for dealing with unmapped reads
        jobIndex.setMaxMemory("5000");
        jobIndex.setQueue(queue);      
        
        return jobIndex;
    }

    public String parameters() {

        String paramCommand = null;
        StringBuilder a = new StringBuilder();

        try {

                if (hasPropertyAndNotNull("additionalStarParams") && !getProperty("additionalStarParams").isEmpty()) {
                    star_aln_params = getProperty("additionalStarParams");
                    a.append(" ");
                    a.append(star_aln_params);
                    a.append(" ");
                }
                
                StringBuilder clipSeq = new StringBuilder();
                if (hasPropertyAndNotNull("r1_adapter_trim") && !getProperty("r1_adapter_trim").isEmpty()) {
                    read1_adapterTrim = getProperty("r1_adapter_trim");
                    clipSeq.append(read1_adapterTrim);
                }
                
                if (hasPropertyAndNotNull("r2_adapter_trim") && !getProperty("r2_adapter_trim").isEmpty()) {
                    read2_adapterTrim = getProperty("r2_adapter_trim");
                    if (clipSeq.length() > 0) {
                        clipSeq.append(" ");
                    }
                       clipSeq.append(read2_adapterTrim);
                }
                
                if (!clipSeq.toString().isEmpty() && clipSeq.toString().length() > 1) {
                    a.append(" ");
                    a.append("--clip3pAdapterSeq ");
                    a.append(clipSeq.toString());
                }
                
                if (this.produceTranscriptomeBam) {
                    a.append(" ");
                    a.append("--quantMode TranscriptomeSAM GeneCounts");
                }
        
                paramCommand = a.toString();
                return paramCommand; 

        } catch (Exception e) {
            System.out.println(e.getMessage());
        }
        return paramCommand;
    }

}
