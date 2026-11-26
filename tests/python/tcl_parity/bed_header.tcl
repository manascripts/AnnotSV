set dir [lindex $argv 0]
source $dir/share/tcl/AnnotSV/AnnotSV-general.tcl
source $dir/share/tcl/AnnotSV/AnnotSV-config.tcl
array set g_AnnotSV [list outputDir [lindex $argv 2] bedFile [lindex $argv 1] svtBEDcol [lindex $argv 3] samplesidBEDcol [lindex $argv 4]]
createBEDinputHeaderFile
addNAinSamplesIDbedCol
puts "bedFile=[file tail $g_AnnotSV(bedFile)] samplesidBEDcol=$g_AnnotSV(samplesidBEDcol)"
