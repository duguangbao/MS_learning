#!perl
#
# Purpose: Demonstrate accessing atoms and their bonding information

use strict;
use warnings;
use MaterialsScript qw(:all);

#A simple utility function to format a Point nicely
sub PrintPoint {
  my ($point) = @_;

  printf "[%2.2f, %2.2f, %2.2f]", $point->X, $point->Y, $point->Z;
}

my $idx = 0;
my $atoms = $Documents{"Phenol.xsd"}->Atoms;

foreach my $atom (@$atoms) {
  print $idx . " Element = " . $atom->ElementName;
  PrintPoint($atom->XYZ);
  print "\n";

  print "\tBond Count = " . $atom->NumBonds . ":";

  if($atom->HasBonds > 0) {
    my $bonds = $atom->Bonds;
    foreach my $bond (@$bonds) {
      print " " . $bond->BondType;
      print " [" . $bond->Atom1->ElementSymbol . "-" .
        $bond->Atom2->ElementSymbol . "]";
    }
  }

  print "\n";
  ++$idx;
}

