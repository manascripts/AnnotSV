set dir [lindex $argv 0]
source $dir/share/tcl/AnnotSV/AnnotSV-general.tcl
source $dir/share/tcl/AnnotSV/AnnotSV-vcf.tcl
array set g_AnnotSV [list outputDir [lindex $argv 2] outputFile x.annotated.tsv SVminSize [lindex $argv 4] bracketedAltMaxSize 10 bracketedAltMode [lindex $argv 5] SVinputInfo [lindex $argv 3] includeCI [lindex $argv 6] missingGTinSamplesid 1 candidateSnvIndelFiles ""]
rename exit _e; set done 0; proc exit {{c 0}} {global g_AnnotSV done; if {!$done} {set f [open $g_AnnotSV(outputDir)/EXIT w]; puts $f $c; close $f}; _e 0}
set bed [VCFsToBED [lindex $argv 1]]
set done 1
file rename $bed $g_AnnotSV(outputDir)/out.bed
set f [open $g_AnnotSV(outputDir)/meta w]
puts $f $VCFheader
foreach k [lsort [array names g_SVLEN]] {puts $f "$k=$g_SVLEN($k)"}
close $f
