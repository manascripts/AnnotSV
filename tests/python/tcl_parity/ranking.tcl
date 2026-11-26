# Runs SVprepareRanking, SVrankingLoss/Gain/INS/INV on every line and achieveSVrankingLoss/Gain,
# then the output order of OrganizeAnnotation (DescendingSortOnElement2 then 1).
# Args: dir input output svinputFile criteriaOverrides(list of key value)
set dir [lindex $argv 0]
foreach m {general ranking} {source $dir/share/tcl/AnnotSV/AnnotSV-$m.tcl}
lassign [lrange $argv 1 end] input output svinput overrides
array set g_AnnotSV [list organism Human SVinputFile $svinput svtTSVcol 5 svtBEDcol 4 ranking 1]
# Default criteria: the "set g_rankingCriteria(...)" lines of AnnotSV-config.tcl
set f [open $dir/share/tcl/AnnotSV/AnnotSV-config.tcl]
foreach l [split [read $f] "\n"] {if {[regexp {^\s*set g_rankingCriteria\([^$]+\)\s+-?[0-9.]+\s*$} $l]} {eval $l}}
close $f
array set g_rankingCriteria $overrides

set f [open $input]
set header [gets $f]
SVprepareRanking $header
set L_AnnotSV_ID {}
while {[gets $f L] >= 0} {
    set Ls [split $L "\t"]
    set AnnotSV_ID [lindex $Ls 0]
    if {[lsearch -exact $L_AnnotSV_ID $AnnotSV_ID] == -1} {lappend L_AnnotSV_ID $AnnotSV_ID}
    if {$g_AnnotSV(ranking)} {
        set SVtype [normalizeSVtype [lindex $Ls 5]]
        if {$SVtype eq "DEL"} {SVrankingLoss $L} elseif {$SVtype eq "DUP"} {SVrankingGain $L
        } elseif {$SVtype eq "INS"} {SVrankingINS $L} elseif {$SVtype eq "INV"} {SVrankingINV $L
        } else {set g_rankingScore($AnnotSV_ID) ""; set g_rankingExplanations($AnnotSV_ID) ""}
    }
}
close $f
set o [open $output w]
puts $o "ranking=$g_AnnotSV(ranking) svtBEDcol=$g_AnnotSV(svtBEDcol) svtTSVcol=$g_AnnotSV(svtTSVcol)"
if {$g_AnnotSV(ranking)} {
    foreach AnnotSV_ID [lsort -unique $L_AnnotSV_ID] {
        set SVtype [normalizeSVtype [lindex [split $AnnotSV_ID "_"] end-1]]
        if {$SVtype eq "DEL"} {achieveSVrankingLoss $AnnotSV_ID} elseif {$SVtype eq "DUP"} {achieveSVrankingGain $AnnotSV_ID}
    }
    set L {}
    foreach AnnotSV_ID $L_AnnotSV_ID {
        if {$g_rankingScore($AnnotSV_ID) eq ""} {set score "-99"} else {set score $g_rankingScore($AnnotSV_ID)}
        lappend L "$AnnotSV_ID $score -1"
    }
    foreach x [lsort -command DescendingSortOnElement1 [lsort -command DescendingSortOnElement2 $L]] {
        set AnnotSV_ID [lindex $x 0]
        set e ""
        if {[info exists g_rankingExplanations($AnnotSV_ID)]} {set e $g_rankingExplanations($AnnotSV_ID)}
        puts $o "$AnnotSV_ID\t$g_rankingScore($AnnotSV_ID)\t$e"
    }
}
close $o
