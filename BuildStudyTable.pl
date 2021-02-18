#!perl
#
# Purpose: Place a molecule at a series of points in a zeolite. Then run
#          a Forcite energy calculation on the structure. The script
#          generates a study table with the input structures and the results.

use strict;
use warnings;
use MaterialsScript qw(:all);

my $doc           = $Documents{"MORCH3Cl.xsd"};

my $newStudyTable = Documents->New("PES.std");

my $chloroMethane = $doc->DisplayRange->Sets("ChloroMethane")->Atoms;

my $calcSheet = $newStudyTable->Sheets->Item(0);

#Setup some Header information in the StudyTable
$calcSheet->ColumnHeading(0) = "Structure";
$calcSheet->ColumnHeading(1) = "Potential Energy";
$calcSheet->ColumnHeading(2) = "Valence Diagonal Energy";

my $forciteEnergy  = Modules->Forcite->Energy;

for(my $zPos = 0; $zPos < 14; ++$zPos) {
  #Move the molecule along the Z axis
  #The initial loop here assumes the molecule is not quite
  #  in the right position.
  $chloroMethane->Translate(Point(X=>0, Y=>0, Z=>-1));

  #For demo only - force the view to update
  $doc->UpdateViews;

  $calcSheet->Cell($zPos, 0) = $doc;

  #Run Forcite
  $forciteEnergy->Run($doc);

  #Populate the study table
  $calcSheet->Cell($zPos, 1) = $doc->PotentialEnergy;
  $calcSheet->Cell($zPos, 2) = $doc->ValenceDiagonalEnergy;
}

$newStudyTable->Save;
