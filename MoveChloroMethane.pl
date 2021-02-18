#!perl
#
# Purpose: Demonstrate moving a molecule through a pore in a zeolite using
#          a named set.

use strict;
use warnings;
use MaterialsScript qw(:all);

my $doc = $Documents{"MORCH3Cl.xsd"};
if (!$doc) {die "no document";}

my $chloroMethane= $doc->DisplayRange->Sets("ChloroMethane")->Atoms;

my $delta = -1;
for(my $x = 0; $x < 4; ++$x) {
    for(my $zPos = 0; $zPos < 14; ++$zPos) {
        #Move the molecule along the Z axis
        #The initial loop here assumes the molecule is not quite
        #  in the right position.
        $chloroMethane->Translate(Point(X=>0, Y=>0, Z=>$delta));

        #For demo only - force the view to update
        $doc->UpdateViews;
    }

    $delta = -$delta;
}
