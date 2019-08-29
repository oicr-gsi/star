package ca.on.oicr.pde.deciders;

import com.google.common.collect.MoreCollectors;
import java.util.Arrays;
import java.util.List;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles;

/**
 *
 * @author mlaszloffy
 */
public class ReadGroupData {

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
        this.tissueType = attrs.stream().map((ca.on.oicr.pde.deciders.FileAttributes a) -> a.getLimsValue(Lims.TISSUE_TYPE)).distinct().collect(MoreCollectors.onlyElement());
        this.lane = attrs.stream().map((a) -> a.getLane().toString()).distinct().collect(MoreCollectors.onlyElement());
        this.RGLB = attrs.stream().map((a) -> a.getLibrarySample()).distinct().collect(MoreCollectors.onlyElement());
        this.RGPU = attrs.stream().map((a) -> a.getSequencerRun() + "_" + this.lane + "_" + a.getBarcode()).distinct().collect(MoreCollectors.onlyElement());
        this.RGSM = attrs.stream().map((ca.on.oicr.pde.deciders.FileAttributes a) -> {
            String rgsm = a.getDonor() + "_" + this.tissueType;
            if (a.getLimsValue(Lims.GROUP_ID) != null && !a.getLimsValue(Lims.GROUP_ID).isEmpty()) {
                rgsm += "_" + a.getLimsValue(Lims.GROUP_ID);
            }
            return rgsm;
        }).distinct().collect(MoreCollectors.onlyElement());
        this.ius_accession = attrs.stream().map((ca.on.oicr.pde.deciders.FileAttributes a) -> a.getOtherAttribute(FindAllTheFiles.Header.IUS_SWA)).distinct().collect(MoreCollectors.onlyElement());
        this.sequencer_run_name = attrs.stream().map((a) -> a.getSequencerRun()).distinct().collect(MoreCollectors.onlyElement());
        this.barcode = attrs.stream().map((a) -> a.getBarcode()).distinct().collect(MoreCollectors.onlyElement());
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
