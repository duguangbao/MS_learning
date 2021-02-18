#!perl

use strict;
use MaterialsScript qw(:all);

# Author: Reinier Akkermans (Accelrys)
# Date: 10 December 2009
# Required: Materials Studio 5.0 or later.

# Converts the input trajectory with atoms and bonds to 
# an output trajectory with beads and connectors.
# Bead types are created with names A, B, C, ...
# The associated structures can be looked up in the bead typing studytable.

my $xtd = $Documents{"Polydimeth_siloxane.xtd"};
my $xtdCG = CoarseGrainTrajectory($xtd);

sub CoarseGrainTrajectory
{
	my ($xtd) = @_;
	my $trj = $xtd->Trajectory;
	
	my $xtdCG = Documents->New($xtd->Name . " CG.xtd");
	my $trjCG = $xtdCG->Trajectory;
	
	my $cg = Tools->CoarseGrainer;
	my $typingDoc = Documents->New("Bead Typing.std");
	
	for (my $frame = 1; $frame <= $trj->NumFrames; ++$frame){
	    $trj->CurrentFrame = $frame;
	    if($frame == 1){
			$cg->CoarseGrain->UpdateTypingDocument($xtd, $typingDoc);
		 }
		 # rename A,B,C, ...
		 for(my $row = 0; $row < $typingDoc->RowCount; $row++){
		 	$typingDoc->Cell($row,1) = chr(65+$row);
		 }
		 my $frameCG = $cg->CoarseGrain->Build($xtd, $typingDoc);
		 $trjCG->AppendFramesFrom($frameCG->Structure);
		 $frameCG->Structure->Delete;
	}
	return $xtdCG;
}