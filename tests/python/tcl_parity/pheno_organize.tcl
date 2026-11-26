# Runs OrganizeAnnotation with the ranking and the PhenoGenius / Exomiser columns. The tool
# results are mocked: g_PhenoGeniusCli and g_Exomiser are read from a TSV file (key, then the
# tab-separated value) as runPhenoGeniusCli / runExomiser would have filled them.
# Same arguments as ranking_organize.tcl, then: mock file, PhenoGeniusCli (0/1), Exomiser (0/1).
# A Tcl error writes "1" in EXIT.
set dir [lindex $argv 0]
foreach m {general ranking config vcf genes write closestgenes omim pathogenicsv pathogenicsnvindel benignsv extann ncbi exomiser phenogenius} {
    source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl
}
set auto_path [linsert $::auto_path 0 $dir/share/tcl/AnnotSV/tcllib]
lassign [lrange $argv 1 end] annDir outDir input build mode cols metrics bcftools extFile rankFiltering svt mockFile pgOn exoOn
array set g_AnnotSV [list annotationsDir $annDir organism Human genomeBuild $build outputDir $outDir \
    outputFile x.annotated.tsv SVinputFile $input svtBEDcol $svt samplesidBEDcol -1 txFile "" \
    annotationMode $mode tx RefSeq SVinputInfo 1 bedtools bedtools bcftools $bcftools SVminSize 50 \
    bracketedAltMaxSize 10 benignAF 0.01 minTotalNumber 500 \
    bracketedAltMode ALL includeCI 1 missingGTinSamplesid 1 candidateSnvIndelFiles "" snvIndelFiles "" \
    cytoband 0 tadAnn 0 gcContentAnn 0 repeatAnn 0 gapAnn 0 segdupAnn 0 ENCODEblacklistAnn 0 \
    geneBasedAnn 1 PhenoGeniusCli $pgOn Exomiser $exoOn ranking 1 candidateGenesFiltering 0 metrics $metrics \
    rankFiltering $rankFiltering promAnn 0 EAann 0 GHann 0 miRNAann 0 ABCann 0 MPRAann 0 overlap 100 \
    extann [list $extFile] outputColHeader [concat AnnotSV_ID Samples_ID $cols] svtTSVcol [expr {$svt eq -1 ? -1 : $svt + 2}] samplesidTSVcol -1 vcf 0]
foreach p {regulatoryElementsAnnotation} {proc $p args {return ""}}
set f [open $mockFile]
foreach l [split [read -nonewline $f] "\n"] {
    set ls [split $l "\t"]
    set v [join [lrange $ls 2 end] "\t"]
    if {[lindex $ls 0] eq "PG"} {set g_PhenoGeniusCli([lindex $ls 1]) $v} else {set g_Exomiser([lindex $ls 1]) $v}
}
close $f
set f [open $dir/share/tcl/AnnotSV/AnnotSV-config.tcl]
foreach l [split [read $f] "\n"] {if {[regexp {^\s*set g_rankingCriteria\([^$]+\)\s+-?[0-9.]+\s*$} $l]} {eval $l}}
close $f
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
