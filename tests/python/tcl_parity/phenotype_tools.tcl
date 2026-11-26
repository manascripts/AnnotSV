# Runs the external-tool procs (Exomiser, PhenoGenius, variantconvert) with mocked tools on PATH,
# then dumps their state. Usage: tclsh phenotype_tools.tcl ROOT WHAT DUMPFILE key value ...
# WHAT: exomiser | phenogenius | variantconvert
set dir [lindex $argv 0]
foreach m {general ncbi exomiser phenogenius variantconvert} {
    source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl
}
set auto_path [linsert $::auto_path 0 $dir/share/tcl/AnnotSV/tcllib]
lassign [lrange $argv 1 2] what dump
array set g_AnnotSV [lrange $argv 3 end]
set out {}
switch $what {
    exomiser {
        checkExomiserInstallation
        lappend out "Exomiser=$g_AnnotSV(Exomiser)" "hpoVersion=[expr {[info exists hpoVersion] ? $hpoVersion : {}}]"
        if {$g_AnnotSV(Exomiser)} {runExomiser $g_AnnotSV(genes) $g_AnnotSV(hpo)}
        lappend out "Exomiser=$g_AnnotSV(Exomiser)"
        foreach k [lsort [array names g_Exomiser]] {lappend out "$k\t$g_Exomiser($k)"}
        foreach g [concat $g_AnnotSV(genes) NOPE] {
            lappend out "$g\t[ExomiserAnnotation $g all]\t[ExomiserAnnotation $g score]"
        }
    }
    phenogenius {
        checkPhenoGeniusCli
        lappend out "PhenoGeniusCli=$g_AnnotSV(PhenoGeniusCli)"
        if {$g_AnnotSV(PhenoGeniusCli)} {
            set L_NCBI_ID {}
            foreach g $g_AnnotSV(genes) {
                set NCBI_ID [searchforNCBIGeneID $g]
                if {$NCBI_ID ne ""} {lappend L_NCBI_ID $NCBI_ID}
            }
            runPhenoGeniusCli $g_AnnotSV(genes) [lsort -unique $L_NCBI_ID] $g_AnnotSV(hpo)
        }
        lappend out "PhenoGeniusCli=$g_AnnotSV(PhenoGeniusCli)"
        foreach k [lsort [array names g_PhenoGeniusCli]] {lappend out "$k\t$g_PhenoGeniusCli($k)"}
        foreach g [concat $g_AnnotSV(genes) NOPE] {
            lappend out "$g\t[PhenoGeniusCliAnnotation $g all]\t[PhenoGeniusCliAnnotation $g specificity]"
        }
    }
    variantconvert {
        checkVariantconvertConfigfile
        runVariantconvert $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile)
    }
}
set f [open $dump w]; puts $f [join $out "\n"]; close $f
