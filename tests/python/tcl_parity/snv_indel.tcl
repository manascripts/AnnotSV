# Runs genesAnnotation and OrganizeAnnotation with the SNV/indel counts (VCFannotation) and the
# compound htz (filteredVCFannotation) columns. A Tcl error writes "1" in EXIT.
set dir [lindex $argv 0]
foreach m {general config vcf filteredVCF genes write closestgenes} {
    source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl
}
set auto_path [linsert $::auto_path 0 $dir/share/tcl/AnnotSV/tcllib]
lassign [lrange $argv 1 end] annDir outDir input mode svt snvFiles snvSamples candFiles candSamples pass
array set g_AnnotSV [list annotationsDir $annDir organism Human genomeBuild GRCh38 outputDir $outDir \
    outputFile x.annotated.tsv SVinputFile $input svtBEDcol $svt samplesidBEDcol -1 txFile "" \
    annotationMode $mode tx RefSeq SVinputInfo 1 bedtools bedtools bcftools bcftools SVminSize 50 \
    bracketedAltMaxSize 10 bracketedAltMode ALL includeCI 1 missingGTinSamplesid 1 \
    candidateSnvIndelFiles $candFiles candidateSnvIndelSamples $candSamples \
    snvIndelFiles $snvFiles snvIndelSamples $snvSamples snvIndelPASS $pass \
    cytoband 0 tadAnn 0 gcContentAnn 0 repeatAnn 0 gapAnn 0 segdupAnn 0 ENCODEblacklistAnn 0 \
    geneBasedAnn 0 PhenoGeniusCli 0 Exomiser 0 ranking 0 candidateGenesFiltering 0 metrics us \
    rankFiltering {1 2 3 4 5 NA} promAnn 0 EAann 0 GHann 0 miRNAann 0 ABCann 0 MPRAann 0 overlap 100 \
    extann "" outputColHeader {AnnotSV_ID Samples_ID} svtTSVcol -1 samplesidTSVcol -1 vcf 0]
foreach p {regulatoryElementsAnnotation SVprepareRanking} {proc $p args {return ""}}
# Not ported in this module's test: always written for Human. The test removes the "@STUB@" fields.
foreach p {pathogenicSVannotation poPathogenicSVannotation pathoSNVindelAnnotation benignSVannotation
           poBenignSVannotation} {proc $p args {return "@STUB@"}}
rename exit _e; set done 0
proc exit {{c 0}} {global g_AnnotSV done; if {$done} {_e 0}; set f [open $g_AnnotSV(outputDir)/EXIT w]; puts $f $c; close $f; _e 0}

if {[catch {
    # configureAnnotSV: compound htz parsing
    if {$g_AnnotSV(candidateSnvIndelFiles) ne ""} {filteredVCFannotation "FULL" "" "" "" ""}
    if {[regexp -nocase "\\.vcf(.gz)?$" $input]} {
        set g_AnnotSV(bedFile) [VCFsToBED $input]
    } else {
        set g_AnnotSV(bedFile) $input
        createBEDinputHeaderFile
        addNAinSamplesIDbedCol
    }
    checkGenesRefSeqFile
    genesAnnotation
    OrganizeAnnotation
} err]} {
    puts stderr $::errorInfo
    set f [open $outDir/EXIT w]; puts $f 1; close $f
}
set done 1
