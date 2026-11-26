# Runs checkGenes*File, genesAnnotation and OrganizeAnnotation with every
# other annotation switched off or stubbed (core output columns only).
set dir [lindex $argv 0]
foreach m {general config vcf genes write} {source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl}
lassign [lrange $argv 1 end] annDir outDir input svt sid txFile mode tx info
array set g_AnnotSV [list annotationsDir $annDir organism Human genomeBuild GRCh38 outputDir $outDir \
    outputFile x.annotated.tsv SVinputFile $input svtBEDcol $svt samplesidBEDcol $sid txFile $txFile \
    annotationMode $mode tx $tx SVinputInfo $info bedtools bedtools SVminSize 50 bracketedAltMaxSize 10 \
    bracketedAltMode ALL includeCI 1 missingGTinSamplesid 1 candidateSnvIndelFiles "" snvIndelFiles "" \
    cytoband 0 tadAnn 0 gcContentAnn 0 repeatAnn 0 gapAnn 0 segdupAnn 0 ENCODEblacklistAnn 0 geneBasedAnn 0 \
    PhenoGeniusCli 0 Exomiser 0 ranking 0 candidateGenesFiltering 0 metrics us rankFiltering {1 2 3 4 5 NA} \
    promAnn 0 EAann 0 GHann 0 miRNAann 0 ABCann 0 MPRAann 0 overlap 100 outputColHeader {AnnotSV_ID Samples_ID} \
    svtTSVcol -1 samplesidTSVcol -1 vcf 0]
foreach p {regulatoryElementsAnnotation SVprepareRanking pathogenicSVannotation poPathogenicSVannotation
           pathoSNVindelAnnotation benignSVannotation poBenignSVannotation} {proc $p args {return ""}}
proc closestGenesAnnotation args {return {}}
rename exit _e; set done 0
proc exit {{c 0}} {global g_AnnotSV done; if {$done} {_e 0}; set f [open $g_AnnotSV(outputDir)/EXIT w]; puts $f $c; close $f; _e 0}

if {[regexp -nocase "\\.vcf(.gz)?$" $input]} {
    set g_AnnotSV(bedFile) [VCFsToBED $input]
} else {
    set g_AnnotSV(bedFile) $input
    createBEDinputHeaderFile
    addNAinSamplesIDbedCol
}
if {$tx eq "RefSeq"} {checkGenesRefSeqFile} else {checkGenesENSEMBLfile}
genesAnnotation
file copy $outDir/x.annotated.tsv.tmp $outDir/fullsplit.txt
set f [open $outDir/genes.txt w]
foreach k [lsort [array names g_Lgenes]] {puts $f "$k=$g_Lgenes($k)"}
close $f
OrganizeAnnotation
set done 1
