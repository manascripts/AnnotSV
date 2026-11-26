# Runs the regulatory elements procs (check*files, regulatoryElementsAnnotation) between
# genesAnnotation and OrganizeAnnotation, with the other annotations switched off or stubbed.
set dir [lindex $argv 0]
foreach m {general config vcf genes write closestgenes regulatoryelements omim candidategenes} {source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl}
lassign [lrange $argv 1 end] annDir outDir input build tx sel1 sel2 report cand
array set g_AnnotSV [list annotationsDir $annDir organism Human genomeBuild $build outputDir $outDir \
    outputFile x.annotated.tsv SVinputFile $input svtBEDcol 4 samplesidBEDcol -1 txFile "" \
    annotationMode both tx $tx SVinputInfo 1 bedtools bedtools SVminSize 50 bracketedAltMaxSize 10 \
    bracketedAltMode ALL includeCI 1 missingGTinSamplesid 1 candidateSnvIndelFiles "" snvIndelFiles "" \
    cytoband 0 tadAnn 0 gcContentAnn 0 repeatAnn 0 gapAnn 0 segdupAnn 0 ENCODEblacklistAnn 0 geneBasedAnn 0 \
    PhenoGeniusCli 0 Exomiser 0 ranking 0 candidateGenesFiltering 0 metrics us rankFiltering {1 2 3 4 5 NA} \
    promoterSize 500 REselect1 $sel1 REselect2 $sel2 REreport $report candidateGenesFile $cand hpo "" \
    overlap 100 outputColHeader {AnnotSV_ID Samples_ID} svtTSVcol -1 samplesidTSVcol -1 vcf 0]
foreach p {SVprepareRanking pathogenicSVannotation poPathogenicSVannotation
           pathoSNVindelAnnotation benignSVannotation poBenignSVannotation} {proc $p args {return ""}}
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
checkGenesENSEMBLfile
foreach p {checkPromoterFile checkEAfiles checkGHfiles checkMiRTargetLinkFiles checkABCfiles checkMPRAfiles} {
    $p
    set f [open $outDir/flags.txt a]; set l {}; foreach k {ABCann EAann GHann MPRAann miRNAann promAnn} {if {[info exists g_AnnotSV($k)]} {lappend l "$k $g_AnnotSV($k)"}}; puts $f "$p [join $l]"; close $f
}
genesAnnotation
set f [open $outDir/hits.txt w]
foreach k [lsort [array names g_HI]] {puts $f "HI $k=$g_HI($k)"}
foreach k [lsort [array names g_TS]] {puts $f "TS $k=$g_TS($k)"}
foreach k [lsort [array names g_re]] {puts $f "RE $k=$g_re($k)"}
close $f
OrganizeAnnotation
set done 1
