
#### Design to find residue pair betweens two proteins in contact


### spike protein chain
mol load pdb 3sci.pdb
set s E  
### receptor protein chain 
set r A  
set sel [atomselect top "( same residue as protein and within 3.5 of chain $s) and not chain $s and  name CA"]
## list of  residues of chain r within 3.5 of chain s
set resname [$sel get resname] 
set out [open "Pair_list.txt" w]
puts $out "CHAIN:$r \t CHAIN:$s"
foreach name $resname {
	set sel2 [atomselect top "(same residue as protein and within 3.5 of (chain $r and resname $name)) and not chain $r and  name CA"]
	## list of  residue of chain s
	set resname2 [ $sel2 get resname] 
	foreach name2 $resname2 { 
		puts "$name \t $name2"
		puts $out "$name \t $name2"

}}
close $out