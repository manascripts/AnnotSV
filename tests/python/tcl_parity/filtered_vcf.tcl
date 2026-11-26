set dir [lindex $argv 0]
source $dir/share/tcl/AnnotSV/AnnotSV-vcf.tcl
source $dir/share/tcl/AnnotSV/AnnotSV-filteredVCF.tcl
array set g_AnnotSV [list SVinputFile x.vcf candidateSnvIndelFiles [lindex $argv 1] candidateSnvIndelSamples [lindex $argv 2] SVminSize 50 snvIndelPASS [lindex $argv 3] formatColNumber 1]
puts "FULL=[filteredVCFannotation FULL {} {} {} {}]|"
foreach k [lsort [array names L_htzPos]] {puts "$k=$L_htzPos($k)"}
puts "samples=$g_AnnotSV(candidateSnvIndelSamples)"
# split line: Line = {x FORMAT s1 s2 ...} + 13 trailing cols
set hdr "A\tFORMAT\t[join $g_AnnotSV(candidateSnvIndelSamples) \t]\tAnnotation_mode"
set line [concat x GT [lrepeat [llength $g_AnnotSV(candidateSnvIndelSamples)] 0/1] [lrepeat 13 z]]
foreach {c s e} {1 1 300000000 1 10000 20000000 2 1 100000000 X 1 200000000} {puts "$c:$s-$e=[filteredVCFannotation $c $s $e $line $hdr]|"}
