# Runs the gene-based check procs (GenCC, NCBI/HGNC, gene aliases, OMIM, morbid, HI, ClinGen,
# LOEUF, ExAC), the extann file list of bin/AnnotSV, genesAnnotation and OrganizeAnnotation
# (gene-based columns, candidate genes filtering) with the other annotations switched off.
set dir [lindex $argv 0]
foreach m {general config vcf genes write closestgenes regulatoryelements omim candidategenes ncbi gencc
           haploinsufficiency clingen loeuf-pLI exac extann} {source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl}
set auto_path [linsert $::auto_path 0 $dir/share/tcl/AnnotSV/tcllib]
lassign [lrange $argv 1 end] annDir outDir input build tx mode metrics candFiltering cand extFile cols
array set g_AnnotSV [list annotationsDir $annDir organism Human genomeBuild $build outputDir $outDir \
    outputFile x.annotated.tsv SVinputFile $input svtBEDcol 4 samplesidBEDcol -1 txFile "" \
    annotationMode $mode tx $tx SVinputInfo 1 bedtools bedtools SVminSize 50 bracketedAltMaxSize 10 \
    bracketedAltMode ALL includeCI 1 missingGTinSamplesid 1 candidateSnvIndelFiles "" snvIndelFiles "" \
    cytoband 0 tadAnn 0 gcContentAnn 0 repeatAnn 0 gapAnn 0 segdupAnn 0 ENCODEblacklistAnn 0 \
    PhenoGeniusCli 0 Exomiser 0 ranking 0 candidateGenesFiltering $candFiltering metrics $metrics \
    rankFiltering {1 2 3 4 5 NA} promAnn 0 EAann 0 GHann 0 miRNAann 0 ABCann 0 MPRAann 0 overlap 100 \
    candidateGenesFile $cand hpo "HP:1" extann "" externalGeneFiles $extFile \
    outputColHeader [concat AnnotSV_ID Samples_ID $cols] svtTSVcol -1 samplesidTSVcol -1 vcf 0]
foreach p {regulatoryElementsAnnotation SVprepareRanking} {proc $p args {return ""}}
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
checkGenesENSEMBLfile
checkGenCCgeneFile
checkNCBIandHGNC
memorizeGeneNameAlias
set f [open $outDir/alias.txt w]
foreach k [lsort [array names g_alias]] {puts $f "$k=$g_alias($k)"}
close $f
checkOMIMfile
checkMorbidfile
checkHIfile
checkClinGenFile
checkLOEUFfile
checkGeneIntoleranceFile
checkCNVintoleranceFile

# bin/AnnotSV
set geneBasedDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based"
foreach annotFile [glob -nocomplain $geneBasedDir/*/*.tsv] {
    if {[regexp "_DGV_samplesInStudies.tsv$" $annotFile]} {continue}
    lappend g_AnnotSV(extann) $annotFile
}
foreach annotFile [glob -nocomplain $geneBasedDir/*/*.tsv.gz] {lappend g_AnnotSV(extann) $annotFile}
foreach annotFile $g_AnnotSV(externalGeneFiles) {lappend g_AnnotSV(extann) $annotFile}
if {$g_AnnotSV(extann) eq ""} {set g_AnnotSV(geneBasedAnn) 0} else {set g_AnnotSV(geneBasedAnn) 1}

genesAnnotation
set f [open $outDir/probe.txt w]
puts $f "hpo=$g_AnnotSV(hpo)"
foreach g [lsort [array names g_Lgenes]] {
    foreach gene [split $g_Lgenes($g) ";"] {
        puts $f "$gene morbid=[isMorbid $gene] candidate=[isCandidate $gene] pheno=[fromOMIMgeneToPhenotype $gene]"
    }
}
foreach id {100000 100001 100002 100003 100004 100005 999} {puts $f "$id [fromOMIMtoPhenotype $id]"}
close $f
OrganizeAnnotation
set done 1
