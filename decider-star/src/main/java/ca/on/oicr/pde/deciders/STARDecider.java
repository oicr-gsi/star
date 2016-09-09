package ca.on.oicr.pde.deciders;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;

/**
 *
 * @author rtahir
 */
public class STARDecider extends OicrDecider {
    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;

    private final String [][] readMateFlags = {{"_R1_","1_sequence.txt",".1.fastq"},{"_R2_","2_sequence.txt",".2.fastq"}};    
    
    private String index_dir;
    private String output_prefix = "./";
    private String output_dir = "seqware-results";
    private String manual_output = "false";
    //private String RGID = "These";
    private String RGLB = "ARE";
    private String RGPL = "Test";
    private String RGPU = "Values";
    private String RGSM = "ADJUST";
    private String RGCM = "";
    private String additionalStarParams = "";
    private String read1_adapterTrim = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG";
    private String read2_adapterTrim = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
    private String numOfThreads = "6";
    private String starMemory   = "16000";
    private String ius_accession;
    private String barcode;
    private String sequencer_run_name;
    private String lane;
    private String queue = "";
    
    private static final String FASTQ_GZ_METATYPE = "chemical/seq-na-fastq-gzip";
    private static final String OICR = "OICR";
  

    public STARDecider() {
        super();
        fileSwaToSmall = new HashMap<String, BeSmall>();
        parser.accepts("ini-file", "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("index-dir", "reference index dir").withRequiredArg();
        parser.accepts("verbose", "Optional: output all SeqWare info.").withRequiredArg();
        parser.accepts("output-prefix", "Optional: the path where the files should be copied to after analysis. output-prefix in INI file.").withRequiredArg();
        parser.accepts("output-dir", "Optional: the folder to put the output into relative to the output-path. Corresponds to output-dir in INI file.").withRequiredArg();
        parser.accepts("manual-output", "Optional*. Set the manual output either to true or false").withRequiredArg();
        parser.accepts("queue", "Queue on SGE cluster.").withRequiredArg();
        
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
        parser.accepts("template-type", "Optional: limit the run to only specified template type").withRequiredArg();
        //Trimming
        parser.accepts("read1-adapter-trim", "Optional: Barcode, default is empty string.").withRequiredArg();
        parser.accepts("read2-adapter-trim", "Optional: Sequencing platform, will be set to production if no value passed.").withRequiredArg();

    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        this.setMetaType(Arrays.asList("chemical/seq-na-fastq-gzip"));
        this.setGroupingStrategy(FindAllTheFiles.Header.FILE_SWA);

        //allows anything defined on the command line to override the defaults here.
        if (this.options.has("index-dir")){
            this.index_dir = options.valueOf("index-dir").toString();
        } else {
            Log.error("index-dir needs to be set");
            System.exit(1);
        }
        if (this.options.has("output-path")) {
            output_prefix = options.valueOf("output-path").toString();
            if (!output_prefix.endsWith("/")) {
                output_prefix += "/";
            }
        }
        
        if (this.options.has("output-dir")) {
            output_dir = options.valueOf("output-dir").toString();
        }

        if (this.options.has("verbose")) {
            Log.setVerbose(true);
        }

        if (this.options.has("manual-output")) {
            this.manual_output = options.valueOf("manual-output").toString();
        }
        
        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        //RG parameters
        if (this.options.has("rg-library")) {
            this.RGLB = options.valueOf("rg-library").toString();
        }
        if (this.options.has("rg-platform")) {
            this.RGPL = options.valueOf("rg-platform").toString();
        }
        if (this.options.has("rg-platform_unit")) {
            this.RGPU = options.valueOf("rg-platform-unit").toString();
            //this.RGID = this.RGPU;
        }
        if (this.options.has("rg-sample-name")) {
            this.RGSM = options.valueOf("rg-sample-name").toString();
        }
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

        ReturnValue val = super.init();

        return val;
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;

            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                try {
                    if (currentRV.getFiles().get(f).getMetaType().equals(FASTQ_GZ_METATYPE)) {
                        metatypeOK = true;
                    }
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                }
            }
            if (!metatypeOK) {
                continue; // Go to the next value
            }

            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);
            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + meta-type
            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();

            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                iusDeetsToRV.put(fileDeets, currentRV);
            } //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(FindAllTheFiles.Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    iusDeetsToRV.put(fileDeets, currentRV);
                }
            }
        }

        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(FindAllTheFiles.Header.FILE_SWA.getTitle())).getGroupByAttribute();
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }
    
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        boolean haveFirstMate = false;
        boolean haveSecondMate = false;

        for (String p : filePaths) {
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }
                if (!haveFirstMate) {haveFirstMate = p.contains("R1");}
                if (!haveSecondMate){haveSecondMate= p.contains("R2");}
                    
                }

            }
        
        if (!haveFirstMate || !haveSecondMate) {
            Log.error("The Decider was not able to find both R1 and R2 fastq files for paired sequencing alignment, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }
        return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }
    
    @Override
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);

        if (this.options.has("template-type")) {
            if (!returnValue.getAttribute(FindAllTheFiles.Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type").equals(this.options.valueOf("template-type"))) {
                return false;
            }
        }   

        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        Log.debug("INI FILE:" + commaSeparatedFilePaths);

        //reset test mode
        if (!this.options.has("test")) {
            this.setTest(false);
        }
        String[] filePaths = commaSeparatedFilePaths.split(",");
        int[] indexes = {0, 1};

        Set fqInputs_end1 = new HashSet();
        Set fqInputs_end2 = new HashSet();
        Set[] fqInputFiles = {fqInputs_end1, fqInputs_end2};
        String fastq_inputs_end_1 = "";
        String fastq_inputs_end_2 = "";

        for (String p : filePaths) {
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }

                for (int i : indexes) {
                    for (int j = 0; j < this.readMateFlags[i].length; j++) {
                        if (p.contains(this.readMateFlags[i][j])) {
                            fqInputFiles[i].add(p);
                            break;
                        }
                    }
                }
            }
        }

        // Format input strings
        if (fqInputFiles[0].size() == 0 || fqInputFiles[1].size() == 0) {
            Log.error("Was not able to retrieve fastq files for either one or two subsets of paired reads, setting mode to test");
            this.setTest(true);
        } else {
            fastq_inputs_end_1 = _join(",", fqInputFiles[0]);
            fastq_inputs_end_2 = _join(",", fqInputFiles[1]);
        }

        Map<String, String> iniFileMap = new TreeMap<String, String>();
        iniFileMap.put("input_file_1", fastq_inputs_end_1);
        iniFileMap.put("input_file_2", fastq_inputs_end_2);
        iniFileMap.put("index_dir", this.index_dir);
        iniFileMap.put("output_prefix", this.output_prefix);
        iniFileMap.put("output_dir", this.output_dir);
        iniFileMap.put("manual_output", this.manual_output);
        //For RG setting
        iniFileMap.put("rg_library", this.RGLB);
        iniFileMap.put("rg_platform", this.RGPL);
        iniFileMap.put("rg_platform_unit", this.RGPU);
        iniFileMap.put("rg_sample_name", this.RGSM);
        iniFileMap.put("rg_organization", this.RGCM);
        iniFileMap.put("additionalStarParams", this.additionalStarParams);

        iniFileMap.put("r1_adapter_trim", this.read1_adapterTrim);
        iniFileMap.put("r2_adapter_trim", this.read2_adapterTrim);
        iniFileMap.put("star_aln_threads", this.numOfThreads);
        iniFileMap.put("star_aln_mem_mb", this.starMemory);

        iniFileMap.put("ius_accession", this.ius_accession);
        iniFileMap.put("sequencer_run_name", this.sequencer_run_name);
        iniFileMap.put("lane", this.lane);
        iniFileMap.put("barcode", this.barcode);

        if (!this.queue.isEmpty()) {
            iniFileMap.put("queue", this.queue);
        } else {
            iniFileMap.put("queue", " ");
        }
        return iniFileMap;
    }

   //Join function
   public static String _join(String separator, Set items) {
       StringBuilder result = new StringBuilder();
       Iterator myItems = items.iterator();
       while(myItems.hasNext()) {
          if (result.length() > 0)
              result.append(separator);

          result.append(myItems.next().toString());
       }

    return result.toString();
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

    private class BeSmall {

        private Date date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String tissueType = null;
        private String path = null;
        private String tubeID = null;
        private String groupID = null;
        private String groupDescription = null;

        public BeSmall(ReturnValue rv) {
            try {
                this.date = format.parse(rv.getAttribute(FindAllTheFiles.Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            this.tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            this.tubeID = fa.getLimsValue(Lims.TUBE_ID);
            if (null == this.tubeID || this.tubeID.isEmpty()) {
                this.tubeID = "NA";
            }
            this.groupID = fa.getLimsValue(Lims.GROUP_ID);
            if (null == this.groupID || this.groupID.isEmpty()) {
                this.groupID = "NA";
            }
            this.groupDescription = fa.getLimsValue(Lims.GROUP_DESC);
            if (null == this.groupDescription || this.groupDescription.isEmpty()) {
                this.groupDescription = "NA";
            }  
            
            lane = fa.getLane().toString();
            RGLB = fa.getLibrarySample();
            RGPU = fa.getSequencerRun() + "_" + lane + "_" + fa.getBarcode();
            RGSM = fa.getDonor() + "_" + tissueType;
            if (!this.groupID.equals("NA")) {
                RGSM = RGSM + "_" + this.groupID;
            }
            
            this.iusDetails = RGLB + RGPU + rv.getAttribute(FindAllTheFiles.Header.FILE_SWA.getTitle());
            ius_accession = rv.getAttribute(FindAllTheFiles.Header.IUS_SWA.getTitle());
            sequencer_run_name = fa.getSequencerRun();
            barcode = fa.getBarcode();
            
            StringBuilder gba = new StringBuilder(fa.getDonor());
            gba.append(":").append(fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE));
            gba.append(":").append(ius_accession);

            String trs = fa.getLimsValue(Lims.TARGETED_RESEQUENCING);
            if (null != trs && !trs.isEmpty()) {
                gba.append(":").append(trs);
            }

            this.groupByAttribute = gba.toString();
            this.path = rv.getFiles().get(0).getFilePath() + "";
        }

        public Date getDate() {
            return this.date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return this.groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getTissueType() {
            return this.tissueType;
        }

        public String getIusDetails() {
            return this.iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return this.path;
        }

        public String getTubeId() {
            return this.tubeID;
        }

        public String getGroupID() {
            return this.groupID;
        }

        public String getGroupDescription() {
            return this.groupDescription;
        }

        public void setPath(String path) {
            this.path = path;
        }
    }
}
