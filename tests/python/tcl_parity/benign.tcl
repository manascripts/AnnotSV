# Runs checkGenesRefSeqFile, checkPathogenicFiles, checkPathoSNVindelFile, checkBenignFiles,
# checkOverlappedGenesBenignFiles, genesAnnotation and OrganizeAnnotation with the benign SV and
# partially overlapped benign SV columns. A Tcl error writes "1" in EXIT.
set dir [lindex $argv 0]
foreach m {general config vcf genes write closestgenes omim pathogenicsv pathogenicsnvindel benignsv} {
    source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl
}
set auto_path [linsert $::auto_path 0 $dir/share/tcl/AnnotSV/tcllib]
lassign [lrange $argv 1 end] annDir outDir input build mode cols metrics bcftools
array set g_AnnotSV [list annotationsDir $annDir organism Human genomeBuild $build outputDir $outDir \
    outputFile x.annotated.tsv SVinputFile $input svtBEDcol 4 samplesidBEDcol -1 txFile "" \
    annotationMode $mode tx RefSeq SVinputInfo 1 bedtools bedtools bcftools $bcftools SVminSize 50 \
    bracketedAltMaxSize 10 benignAF 0.01 minTotalNumber 500 \
    bracketedAltMode ALL includeCI 1 missingGTinSamplesid 1 candidateSnvIndelFiles "" snvIndelFiles "" \
    cytoband 0 tadAnn 0 gcContentAnn 0 repeatAnn 0 gapAnn 0 segdupAnn 0 ENCODEblacklistAnn 0 \
    geneBasedAnn 0 PhenoGeniusCli 0 Exomiser 0 ranking 0 candidateGenesFiltering 0 metrics $metrics \
    rankFiltering {1 2 3 4 5 NA} promAnn 0 EAann 0 GHann 0 miRNAann 0 ABCann 0 MPRAann 0 overlap 100 \
    extann "" outputColHeader [concat AnnotSV_ID Samples_ID $cols] svtTSVcol -1 samplesidTSVcol -1 vcf 0]
foreach p {regulatoryElementsAnnotation SVprepareRanking} {proc $p args {return ""}}
rename exit _e; set done 0
proc exit {{c 0}} {global g_AnnotSV done; if {$done} {_e 0}; set f [open $g_AnnotSV(outputDir)/EXIT w]; puts $f $c; close $f; _e 0}

if {[catch {
    if {[regexp -nocase "\\.vcf(.gz)?$" $input]} {
        set g_AnnotSV(bedFile) [VCFsToBED $input]
    } else {
        set g_AnnotSV(bedFile) $input
        createBEDinputHeaderFile
        addNAinSamplesIDbedCol
    }
    checkGenesRefSeqFile
    checkPathogenicFiles
    checkPathoSNVindelFile
    checkBenignFiles
    checkOverlappedGenesBenignFiles
    set f [open $outDir/colheader.txt w]; puts $f [join $g_AnnotSV(outputColHeader) "\n"]; close $f
    genesAnnotation
    OrganizeAnnotation
} err]} {
    puts stderr $::errorInfo
    set f [open $outDir/EXIT w]; puts $f 1; close $f
}
set done 1
