# Runs the region check procs (GC content, Repeat, SegDup, Gap, ENCODE blacklist, cytoband, TAD,
# COSMIC, users BED), genesAnnotation and OrganizeAnnotation with the region columns switched on.
set dir [lindex $argv 0]
foreach m {general ranking config vcf genes write closestgenes cytoband gap segdup repeat encodeblacklist tad
           userBED gccontent cosmic} {source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl}
set auto_path [linsert $::auto_path 0 $dir/share/tcl/AnnotSV/tcllib]
lassign [lrange $argv 1 end] annDir outDir input build mode metrics overlap reciprocal cols
array set g_AnnotSV [list annotationsDir $annDir organism Human genomeBuild $build outputDir $outDir \
    outputFile x.annotated.tsv SVinputFile $input svtBEDcol 4 samplesidBEDcol -1 txFile "" \
    annotationMode $mode tx RefSeq SVinputInfo 1 bedtools bedtools SVminSize 50 bracketedAltMaxSize 10 \
    bracketedAltMode ALL includeCI 1 missingGTinSamplesid 1 candidateSnvIndelFiles "" snvIndelFiles "" \
    geneBasedAnn 0 PhenoGeniusCli 0 Exomiser 0 ranking 0 candidateGenesFiltering 0 metrics $metrics \
    rankFiltering {1 2 3 4 5 NA} promAnn 0 EAann 0 GHann 0 miRNAann 0 ABCann 0 MPRAann 0 overlap $overlap \
    reciprocal $reciprocal breakpointProximity 100 extann "" \
    outputColHeader [concat AnnotSV_ID Samples_ID $cols] svtTSVcol -1 samplesidTSVcol -1 vcf 0]
foreach p {regulatoryElementsAnnotation} {proc $p args {return ""}}
# Not ported yet: always written for Human, whatever the header. The test removes the "@STUB@" fields.
foreach p {pathogenicSVannotation poPathogenicSVannotation pathoSNVindelAnnotation benignSVannotation
           poBenignSVannotation} {proc $p args {return "@STUB@"}}
rename exit _e; set done 0
proc exit {{c 0}} {global g_AnnotSV done; if {$done} {_e 0}; set f [open $g_AnnotSV(outputDir)/EXIT w]; puts $f $c; close $f; _e 0}

if {[regexp -nocase "\\.vcf(.gz)?$" $input]} {
    set g_AnnotSV(bedFile) [VCFsToBED $input]
} else {
    set g_AnnotSV(bedFile) $input
    createBEDinputHeaderFile
    addNAinSamplesIDbedCol
}
checkGenesRefSeqFile
# bin/AnnotSV order
checkFASTAfiles
checkRepeatFile
checkSegDupFile
checkGapFile
checkENCODEblacklistFile
checkCytoband
checkTADfiles
checkCOSMICfile
checkUsersBED
genesAnnotation
OrganizeAnnotation
set done 1
