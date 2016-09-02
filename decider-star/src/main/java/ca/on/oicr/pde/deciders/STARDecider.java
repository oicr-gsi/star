package com.github.seqware;

import ca.on.oicr.pde.deciders.FileAttributes;
import ca.on.oicr.pde.deciders.Lims;
import ca.on.oicr.pde.deciders.OicrDecider;
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

    private String [][] readMateFlags = {{"_R1_","1_sequence.txt",".1.fastq"},{"_R2_","2_sequence.txt",".2.fastq"}};    
    
    private String index = "${workflow_bundle_dir}/Workflow_Bundle_${workflow-directory-name}/${version}/data/0.6.2/hg19_random.fa";
    private String run_ends = "2";
    private String output_prefix = "./";
    private String output_dir = "seqware-results";
    private String outputFileName = "";
    private String setManualpath = "";
    private String RGID = "These";
    private String RGLB = "ARE";
    private String RGPL = "Test";
    private String RGPU = "Values";
    private String RGSM = "ADJUST";
    private String additionalPicardParams = "";
    private String adapter_Trimming_activated = "";
    private String read1_adapterTrim = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG";
    private String read2_adapterTrim = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
    private String readTrimming = "";
    private String numOfThreads = "8";
    private String pairingAccuracy = "";
    private String maxInsertSize = "";
    private String readGroup = "";
    private String bwa_aln_params = "";
    private String bwa_sampe_params = "";
    private String ius_accession;
    private Object sequencer_run_name;
    private Object lane;
  

    public STARDecider() {
        super();
        parser.accepts("ini-file", "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("index", "reference ");
        parser.accepts("verbose", "Optional: output all SeqWare info.").withRequiredArg();
        parser.accepts("output-prefix", "Optional: the path where the files should be copied to after analysis. output-prefix in INI file.").withRequiredArg();
        parser.accepts("output-dir", "Optional: the folder to put the output into relative to the output-path. Corresponds to output-dir in INI file.").withRequiredArg();
        parser.accepts("output-filename", "Optional: Template type for grouping samples.").withRequiredArg();
        parser.accepts("set-manual-path", "Optional: colorspace for Novoalign analysis, default 0.").withRequiredArg();
        parser.accepts("run-ends","Run ends will define if it is Single-End(1) or Paired-End(2) experiment, default 2.").withRequiredArg();
        //bwa aln
        parser.accepts("bwa-read-trimming", "Optional: Picard merge slots, default 1.").withRequiredArg();
        parser.accepts("bwa-threads", "Optional: Picard merge slots, default 1.").withRequiredArg();
        parser.accepts("bwa-pairing-accuracy", "Optional: Picard merge slots, default 1.").withRequiredArg();
        parser.accepts("bwa-additional-alignment-parameters", "Optional: Picard merge slots, default 1.").withRequiredArg();
        //bwa sampe
        parser.accepts("bwa-max-insert", "Optional: Picard merge slots, default 1.").withRequiredArg();
        parser.accepts("bwa-read-group", "Optional: Picard merge slots, default 1.").withRequiredArg();
        parser.accepts("bwa-sampe-parameters", "Optional: Picard merge slots, default 1.").withRequiredArg();
        //Picard parameters
        parser.accepts("RGID", "Optional: Picard threads, default 1.").withRequiredArg();
        parser.accepts("RGLB", "Optional: Picard slots, default 1.").withRequiredArg();
        parser.accepts("RGPL", "Optional: Picard memory, default 3000.").withRequiredArg();
        parser.accepts("RGPU", "Optional: Picard merge slots, default 1.").withRequiredArg();
        parser.accepts("RGSM", "Optional: Picard merge slots, default 1.").withRequiredArg();
        parser.accepts("additional-picard-params", "Sample name, default sample.").withRequiredArg();
        //Hotfix addition
        parser.accepts("adapter-trimming");
        parser.accepts("read1-adapter-trim", "Optional: Barcode, default is empty string.").withRequiredArg();
        parser.accepts("read2-adapter-trim", "Optional: Sequencing platform, will be set to production if no value passed.").withRequiredArg();

    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
        //this.setHeader(Header.IUS_SWA);
        this.setMetaType(Arrays.asList("chemical/seq-na-fastq", "chemical/seq-na-fastq-gzip"));

        //allows anything defined on the command line to override the defaults here.
        if (this.options.has("index")){
            this.index = options.valueOf("index").toString();
        }
        if (this.options.has("output-path")) {
            output_prefix = options.valueOf("output-path").toString();
            if (!output_prefix.endsWith("/")) {
                output_prefix += "/";
            }
        }
        if (this.options.has("run-ends")) {
	    String runEnds = options.valueOf("run-ends").toString();
            //verifies that there are two input files
            //verification must take place because STAR decider is not designed for single reads
            if (!runEnds.equals("2")) {
	        Log.error("You passed run-ends parameter " + runEnds + ", but this decider irecognizes only 1 (single reads) and 2 (paired-reads) options");
                System.exit(1);
	    }
            this.run_ends  = options.valueOf("run-ends").toString();
        }
        if (this.options.has("output-dir")) {
            output_dir = options.valueOf("output-dir").toString();
        }
        if (this.options.has("output-filename")) {
            outputFileName = options.valueOf("output-filename").toString();
        }

        if (this.options.has("verbose")) {
            Log.setVerbose(true);
        }

        if (this.options.has("set-manual-path")) {
            this.setManualpath = options.valueOf("set-manual-path").toString();
        }

        //Picard parameters
        if (this.options.has("RGID")) {
            this.RGID = options.valueOf("RGID").toString();
        }
        if (this.options.has("RGLB")) {
            this.RGLB = options.valueOf("RGLB").toString();
        }
        if (this.options.has("RGPL")) {
            this.RGPL = options.valueOf("RGPL").toString();
        }
        if (this.options.has("RGPU")) {
            this.RGPU = options.valueOf("RGPU").toString();
        }
        if (this.options.has("RGSM")) {
            this.RGSM = options.valueOf("RGSM").toString();
        }
        if (this.options.has("additional-picard-params")) {
            this.additionalPicardParams = options.valueOf("additional-picard-params").toString();
        }

        //bwa
        if (this.options.has("bwa-read-trimming")) {
            this.readTrimming = options.valueOf("bwa-read-trimming").toString();
        }
        if (this.options.has("bwa-threads")) {
            this.numOfThreads = options.valueOf("bwa-threads").toString();
        }
        if (this.options.has("bwa-pairing-accuracy")) {
            this.pairingAccuracy = options.valueOf("bwa-pairing-accuracy").toString();
        }
        if (this.options.has("bwa-additional-alignment-parameters")) {
            this.bwa_aln_params = options.valueOf("bwa-additional-alignment-parameters").toString();
        }

        if (this.options.has("bwa-max-insert")) {
            this.maxInsertSize = options.valueOf("bwa-max-insert").toString();
        }
        if (this.options.has("bwa-read-group")) {
            this.readGroup = options.valueOf("bwa-read-group").toString();
        }
        if (this.options.has("bwa-sampe-parameters")) {
            this.bwa_sampe_params = options.valueOf("bwa-sampe-parameters").toString();
        }

        if (this.options.has("adapter-trimming")) {
            this.adapter_Trimming_activated = options.valueOf("adapter-trimming").toString();
        }
        if (this.options.has("read1-adapter-trim")) {
            this.read1_adapterTrim = options.valueOf("read1-adapter-trim").toString();
        }
        if (this.options.has("read2-adapter-trim")) {
            this.read2_adapterTrim = options.valueOf("read2-adapter-trim").toString();
        }

        ReturnValue val = super.init();

        return val;
    }

    protected String handleGroupByAttribute(String attribute, String template, String group_id) {
        //group by parent name, group_id  and template type
        String[] parentNames = attribute.split(":");
        String groupBy = "";
        String[] myFilters = {parentNames[parentNames.length - 1], template, group_id};

        for (int i = 0; i < myFilters.length; i++) {
            if (null != myFilters[i]) {
                if (groupBy.length() > 1) {
                    groupBy = groupBy.concat(":" + myFilters[i]);
                } else {
                    groupBy = groupBy.concat(myFilters[i]);
                }
            }
        }
        return groupBy;
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);

        if (this.options.has("template-type")) {
            if (!returnValue.getAttribute(FindAllTheFiles.Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type").equals(this.options.valueOf("template-type"))) {
                return false;
            }
        }
        
        FileAttributes rv = new  FileAttributes(returnValue, returnValue.getFiles().get(0));
        this.RGID = rv.getLimsValue(Lims.GROUP_ID);
        this.RGLB = rv.getLibrarySample() + rv.getLimsValue(Lims.GROUP_ID);
        this.RGPL = "illumina";
        this.RGPU = rv.getSequencerRun() 
                                + "-" 
                                + rv.getBarcode() 
                                + "_" 
                                + lane;
        this.RGSM= rv.getDonor() + rv.getLimsValue(Lims.GROUP_ID);
//        //Get additional metadata
//        if (null != this.ius_accession) {
//            this.ius_accession = this.ius_accession + "," + returnValue.getAttribute(FindAllTheFiles.Header.IUS_SWA.getTitle());
//        } else {
//            this.ius_accession = returnValue.getAttribute(FindAllTheFiles.Header.IUS_SWA.getTitle());
//        }
//
//        if (null != this.sequencer_run_name) {
//            this.sequencer_run_name = this.sequencer_run_name + "," + returnValue.getAttribute(FindAllTheFiles.Header.SEQUENCER_RUN_NAME.getTitle());
//        } else {
//            this.sequencer_run_name = returnValue.getAttribute(FindAllTheFiles.Header.SEQUENCER_RUN_NAME.getTitle());
//        }
//
//        if (null != this.lane) {
//            this.lane = this.lane + "," + returnValue.getAttribute(FindAllTheFiles.Header.LANE_NUM.getTitle());
//        } else {
//            this.lane = returnValue.getAttribute(FindAllTheFiles.Header.LANE_NUM.getTitle());
//        }

        return super.checkFileDetails(returnValue, fm);
    }

     protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        Log.debug("INI FILE:" + commaSeparatedFilePaths);
        String skipFile = "";
       //reset test mode
        if (!this.options.has("test")) {
            this.setTest(false);
        }
  
        Set fqInputs_end1 = new HashSet();
	Set fqInputs_end2 = new HashSet();
	Set [] fqInputFiles = {fqInputs_end1,fqInputs_end2};
	String fastq_inputs_end_1 = "";
	String fastq_inputs_end_2 = "";

        if (commaSeparatedFilePaths.contains(",")) {
	  String [] fqFilesArray = commaSeparatedFilePaths.split(",");
	  Set fqFilesSet = new HashSet(Arrays.asList(fqFilesArray));
	  Iterator fqFiles = fqFilesSet.iterator();
          int [] indexes = {0,1};

	  while(fqFiles.hasNext()) {
		String file = fqFiles.next().toString();
		for (int i:indexes) {
		   for (int j=0; j < readMateFlags[i].length; j++) {
		     if (file.contains(readMateFlags[i][j])) {
			fqInputFiles[i].add(file);
			break;
		     }
		   }
		}
	 }

	 if (fqInputFiles[0].size() == 0 || fqInputFiles[1].size() == 0) {
	    Log.error("Was not able to retrieve fastq files for either one or two subsets of paired reads, setting mode to test");
	    this.setTest(true);
	 } else {
	    fastq_inputs_end_1 = _join(",",fqInputFiles[0]);
	    fastq_inputs_end_2 = _join(",",fqInputFiles[1]);
	 }
        } else {
          
          fastq_inputs_end_1 = commaSeparatedFilePaths;
          fastq_inputs_end_2 = commaSeparatedFilePaths;
        }

	Map<String, String> iniFileMap = new TreeMap<String, String>();
	iniFileMap.put("input_file_1",fastq_inputs_end_1);
	iniFileMap.put("input_file_2",fastq_inputs_end_2);
	iniFileMap.put("input_reference", this.index);
        iniFileMap.put("output_prefix",this.output_prefix);
	iniFileMap.put("output_dir", this.output_dir);
        iniFileMap.put("outputFileName", this.outputFileName);
        iniFileMap.put("setManualpath", this.setManualpath);
        //For picard
        iniFileMap.put("RGID", this.RGID);
        iniFileMap.put("RGLB", this.RGLB);
        iniFileMap.put("RGPL", this.RGPL);
        iniFileMap.put("RGPU", this.RGPU);
	iniFileMap.put("RGSM", this.RGSM);
	iniFileMap.put("additionalPicardParams", this.additionalPicardParams);
        
        iniFileMap.put("adapter_Trimming_activated", this.adapter_Trimming_activated);
        iniFileMap.put("read1_adapterTrim", this.read1_adapterTrim);
        iniFileMap.put("read2_adapterTrim", this.read2_adapterTrim);
        //For Picard
        iniFileMap.put("readTrimming", this.readTrimming);
        iniFileMap.put("numOfThreads", this.numOfThreads);
        iniFileMap.put("pairingAccuracy", this.pairingAccuracy);
        iniFileMap.put("bwa_aln_params", this.bwa_aln_params);
        
        iniFileMap.put("maxInsertSize", this.maxInsertSize);
        iniFileMap.put("readGroup", this.readGroup);
        iniFileMap.put("bwa_sampe_params", this.bwa_sampe_params);
        //For Read Group
       
	//Hotfix addition
//        iniFileMap.put("queue", this.queue);
//
//	if (this.run_ends.equals("2")) {
//		iniFileMap.put("barcode", this.barcode + "," + this.barcode);
//		iniFileMap.put("ius_accession", _getLastN(this.ius_accession,2));
//		iniFileMap.put("sequencer_run_name", _getLastN(this.sequencer_run_name,2));
//		iniFileMap.put("lane", _getLastN(this.lane,2));
//	} else {
//		iniFileMap.put("barcode", this.barcode);
//		iniFileMap.put("ius_accession", _getLastN(this.ius_accession,1));
//		iniFileMap.put("sequencer_run_name", _getLastN(this.sequencer_run_name,1));
//		iniFileMap.put("lane", _getLastN(this.lane,1));
//	}

        return iniFileMap;
	}

   //Join function
   public static String _join(String separator, Set items) {
       StringBuffer result = new StringBuffer();
       Iterator myItems = items.iterator();
       while(myItems.hasNext()) {
          if (result.length() > 0)
              result.append(separator);

          result.append(myItems.next().toString());
       }

    return result.toString();
    }

    //Element extractor - create array from comma-separated list and return comma-joined last n elements
    public static String _getLastN(String input,int last) {
       String [] elements = input.split(",");
       int start = elements.length - last;
       if (start < 0) {
        Log.error("Attempt to extract more elements than there are in the list " + input + " gets " + elements.length + "Elements");
        return elements[elements.length -1 ]; // return just the last one
       }

       String result = elements[start];
       for (int i = start + 1; i < elements.length; i++) {
	 result = result + "," + elements[i];
	 }

     return result;
    }
    
    public static void main(String args[]){
  
        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(STARDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));
          
    }
}
