#!perl
#
# Purpose: This script builds a 'correct' methane molecule from scratch

use strict;
use warnings;
use MaterialsScript qw(:all);

#Utility function to create an atom of a specific type in a doc
#  Optionally creates a bond to that atom
sub CreateAtom {
    # Bond information ($bondTo, $bondType) is optional
    my ($doc, $element, $pos, $bondTo, $bondType) = @_;
  
    my $newAtom = $doc->CreateAtom($element, $pos);
  
    # Only build the bond if the bond information is present
    if($bondTo && $bondType) {
        $doc->CreateBond($newAtom, $bondTo, $bondType);
    }
  
    return $newAtom;
}

my $doc = Documents->New("Methane.xsd");

#To demonstrate Forcite Minimization put the H atoms in a planar cross
my $centralC = CreateAtom($doc, "C", Point(X => 1, Y => 1, Z => 1));
CreateAtom($doc, "H", Point(X => 1, Y => 0, Z => 1), $centralC, "Single");
CreateAtom($doc, "H", Point(X => 2, Y => 1, Z => 1), $centralC, "Single");
CreateAtom($doc, "H", Point(X => 1, Y => 2, Z => 1), $centralC, "Single");
CreateAtom($doc, "H", Point(X => 0, Y => 1, Z => 1), $centralC, "Single");

#Minimize the structure
Modules->Forcite->GeometryOptimization->Run($doc);

$doc->Save;
