package ca.on.oicr.pde.deciders;

import com.google.common.collect.MoreCollectors;
import com.google.common.collect.Sets;
import java.util.*;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.module.ReturnValue.ExitStatus;
import net.sourceforge.seqware.common.util.Log;

/**
 *
 * @author rtahir
 */
public class STARDecider extends OicrDecider {
    private String index_dir;
    private String produce_transcriptome_bam = "true";
    private String RGCM = "";
    private String additionalStarParams = "";
    private String read1_adapterTrim = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
    private String read2_adapterTrim = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
    private String numOfThreads = "6";
    private String starMemory   = "16000";
    private String queue = "";
    
    private static final String OICR = "OICR";
    private static final String ILLUMINA = "Illumina"; //If we don't have this passed as parameter, we assume Illumina
    private Set<String> allowedTemplateTypes;
    
    private String inputFile1;
    private String inputFile2;
    private ReadGroupData readGroupDataForWorkflowRun;
    
    public STARDecider() {
        super();
        parser.accepts("ini-file", "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("index-dir", "reference index dir").withRequiredArg();
        parser.accepts("produce-transcriptome-bam", "Optional*. Set either to true or false for production of additional files").withRequiredArg();
        
        //star aln
        parser.accepts("star-aln-threads", "Optional: STAR threads, default is 6.").withRequiredArg();
        parser.accepts("star-aln-mem-mb", "Optional: STAR allocated memory Mb, default is 16000.").withRequiredArg();
        parser.accepts("additionalStarParams", "Optional: Star additional parameters").withRequiredArg();
        //RG parameters
        parser.accepts("rg-library",  "Optional: RG Library, (LB).").withRequiredArg();
        parser.accepts("rg-platform", "Optional: RG Platform, default illumina.").withRequiredArg();
        parser.accepts("rg-platform-unit", "Optional: RG Platform unit (PU).").withRequiredArg();
        parser.accepts("rg-sample-name", "Optional: RG Sample name (SM).").withRequiredArg();
        parser.accepts("rg-organization", "Optional: RG Organization (CM).").withRequiredArg();
        parser.accepts("template-type", "Optional: limit the run to only specified template type(s) (comma separated list).").withRequiredArg();
        //Trimming
        parser.accepts("r1-adapter-trim", "Optional: Barcode, default is empty string.").withRequiredArg();
        parser.accepts("r2-adapter-trim", "Optional: Sequencing platform, will be set to production if no value passed.").withRequiredArg();

    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList("chemical/seq-na-fastq-gzip"));
        this.setHeadersToGroupBy(Arrays.asList(FindAllTheFiles.Header.IUS_SWA));

        //allows anything defined on the command line to override the defaults here.
        if (this.options.has("index-dir")){
            this.index_dir = options.valueOf("index-dir").toString();
        } else {
            Log.error("index-dir needs to be set");
            System.exit(1);
        }
        
        if (this.options.has("produce-transcriptome-bam")) {
            this.produce_transcriptome_bam = options.valueOf("produce-transcriptome-bam").toString();
        }

        //RG parameters populated at the end     
        if (this.options.has("rg-organization")) {
            this.RGCM = options.valueOf("rg-organization").toString();
        } else {
            this.RGCM = OICR;
        }
        
        //star
        if (this.options.has("star-threads")) {
            this.numOfThreads = options.valueOf("star-threads").toString();
        }
        if (this.options.has("star-aln-mem-mb")) {
            this.starMemory = options.valueOf("star-aln-mem-mb").toString();
        }
        if (this.options.has("r1-adapter-trim")) {
            this.read1_adapterTrim = options.valueOf("r1-adapter-trim").toString();
        }
        if (this.options.has("r2-adapter-trim")) {
            this.read2_adapterTrim = options.valueOf("r2-adapter-trim").toString();
        }
        if (this.options.has("template-type")) {
            String templateTypeArg = this.options.valueOf("template-type").toString();
            allowedTemplateTypes = Sets.newHashSet(templateTypeArg.split(","));
        }

        ReturnValue val = super.init();

        return val;
    }

    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        this.inputFile1 = null;
        this.inputFile2 = null;

        String[] filePaths = commaSeparatedFilePaths.split(",");
        if (filePaths.length != 2) {
            Log.error("This Decider supports only cases where we have only 2 files per lane, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }

        String[] fqFilesArray = commaSeparatedFilePaths.split(",");
        for (String file : fqFilesArray) {
            int mate = idMate(file);
            switch (mate) {
                case 1:
                    if (this.inputFile1 != null) {
                        Log.error("More than one file found for read 1: " + inputFile1 + ", " + file);
                        return new ReturnValue(ExitStatus.INVALIDFILE);
                    }
                    this.inputFile1 = file;
                    break;
                case 2:
                    if (this.inputFile2 != null) {
                        Log.error("More than one file found for read 2: " + inputFile2 + ", " + file);
                        return new ReturnValue(ExitStatus.INVALIDFILE);
                    }
                    this.inputFile2 = file;
                    break;
                default:
                    Log.error("Cannot identify " + file + " end (read 1 or 2)");
                    return new ReturnValue(ExitStatus.INVALIDFILE);
            }
        }

        if (inputFile1 == null || inputFile2 == null) {
            Log.error("The Decider was not able to find both R1 and R2 fastq files for paired sequencing alignment, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }
        
        readGroupDataForWorkflowRun = new ReadGroupData(files.get(inputFile1), files.get(inputFile2));
        
        return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }
        
    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);

        if (allowedTemplateTypes != null) {
            String currentTemplateType = returnValue.getAttribute(FindAllTheFiles.Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
            if (!allowedTemplateTypes.contains(currentTemplateType)) {
                return false;
            }
        }

        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        Log.debug("INI FILE:" + commaSeparatedFilePaths);

        Map<String, String> iniFileMap = super.modifyIniFile(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        iniFileMap.put("input_file_1", inputFile1);
        iniFileMap.put("input_file_2", inputFile2);
        iniFileMap.put("index_dir", this.index_dir);
        iniFileMap.put("produce_transcriptome_bam", this.produce_transcriptome_bam);
        
        //For RG setting
        String RGLB = this.options.has("rg-library") ? options.valueOf("rg-library").toString() : readGroupDataForWorkflowRun.getRGLB();
        String RGPL = this.options.has("rg-platform") ? options.valueOf("rg-platform").toString() : ILLUMINA; 
        String RGPU = this.options.has("rg-platform_unit") ? options.valueOf("rg-platform-unit").toString() : readGroupDataForWorkflowRun.getRGPU();
        String RGSM = this.options.has("rg-sample-name") ? options.valueOf("rg-sample-name").toString() : readGroupDataForWorkflowRun.getRGSM();

        iniFileMap.put("rg_library", RGLB);
        iniFileMap.put("rg_platform", RGPL);
        iniFileMap.put("rg_platform_unit", RGPU);
        iniFileMap.put("rg_sample_name", RGSM);
        iniFileMap.put("rg_organization", this.RGCM);
        iniFileMap.put("additionalStarParams", this.additionalStarParams);

        iniFileMap.put("r1_adapter_trim", this.read1_adapterTrim);
        iniFileMap.put("r2_adapter_trim", this.read2_adapterTrim);
        iniFileMap.put("star_aln_threads", this.numOfThreads);
        iniFileMap.put("star_aln_mem_mb", this.starMemory);

        iniFileMap.put("ius_accession", readGroupDataForWorkflowRun.getIus_accession());
        iniFileMap.put("sequencer_run_name", readGroupDataForWorkflowRun.getSequencer_run_name());
        iniFileMap.put("lane", readGroupDataForWorkflowRun.getLane());
        iniFileMap.put("barcode", readGroupDataForWorkflowRun.getBarcode());

        return iniFileMap;
    }
   
    public static void main(String args[]) {

        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(STARDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }

    private class ReadGroupData {

        private final String tissueType;
        private final String RGLB;
        private final String RGPU;
        private final String RGSM;
        private final String ius_accession;
        private final String sequencer_run_name;
        private final String barcode;
        private final String lane;

        public ReadGroupData(FileAttributes... inputFileAttrs) {
            List<FileAttributes> attrs = Arrays.asList(inputFileAttrs);
            this.tissueType = attrs.stream().map(a -> a.getLimsValue(Lims.TISSUE_TYPE)).distinct().collect(MoreCollectors.onlyElement());
            this.lane = attrs.stream().map(a -> a.getLane().toString()).distinct().collect(MoreCollectors.onlyElement());
            this.RGLB = attrs.stream().map(a -> a.getLibrarySample()).distinct().collect(MoreCollectors.onlyElement());
            this.RGPU = attrs.stream().map(a -> a.getSequencerRun() + "_" + this.lane + "_" + a.getBarcode()).distinct().collect(MoreCollectors.onlyElement());
            this.RGSM = attrs.stream().map(a -> {
                String rgsm = a.getDonor() + "_" + this.tissueType;
                if (a.getLimsValue(Lims.GROUP_ID) != null && !a.getLimsValue(Lims.GROUP_ID).isEmpty()) {
                    rgsm += "_" + a.getLimsValue(Lims.GROUP_ID);
                }
                return rgsm;
            }).distinct().collect(MoreCollectors.onlyElement());
            this.ius_accession = attrs.stream().map(a -> a.getOtherAttribute(FindAllTheFiles.Header.IUS_SWA)).distinct().collect(MoreCollectors.onlyElement());
            this.sequencer_run_name = attrs.stream().map(a -> a.getSequencerRun()).distinct().collect(MoreCollectors.onlyElement());
            this.barcode = attrs.stream().map(a -> a.getBarcode()).distinct().collect(MoreCollectors.onlyElement());
        }

        public String getTissueType() {
            return this.tissueType;
        }

        public String getRGLB() {
            return RGLB;
        }

        public String getRGPU() {
            return RGPU;
        }

        public String getRGSM() {
            return RGSM;
        }

        public String getIus_accession() {
            return ius_accession;
        }

        public String getSequencer_run_name() {
            return sequencer_run_name;
        }

        public String getBarcode() {
            return barcode;
        }

        public String getLane() {
            return lane;
        }
    }
}
