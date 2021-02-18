#!perl
#
# Purpose: Generate a list of all possible bond lengths/angles in a crystal
#          Any symmetry replicates are ignored

use strict;
use warnings;
use MaterialsScript qw(:all);

use Math::Trig;

#Compute a vector as the delta between point1 & point2
sub Vector {
  my ($point1, $point2) = @_;

  return Point(X => $point1->X - $point2->X, 
               Y => $point1->Y - $point2->Y,
               Z => $point1->Z - $point2->Z);
}

#Calculate the length of a vector
sub Length {
  my ($vector) = @_;

  return sqrt($vector->X * $vector->X + $vector->Y * $vector->Y + $vector->Z * $vector->Z);
}

#Calculate the dot product of two vectors
sub DotProduct {
  my ($vector1, $vector2) = @_;

  return $vector1->X * $vector2->X + $vector1->Y * $vector2->Y + $vector1->Z * $vector2->Z;
}

#Calculate the angle between three points, point1 is at the vertex of the angle
sub CalcAngle {
  my ($point1, $point2, $point3) = @_;

  my $v1 = Vector($point2, $point1);
  my $v2 = Vector($point3, $point1);

  return acos(DotProduct($v1, $v2)/(Length($v1) * Length($v2)));
}

#main

#Define a format for the output
my ($bondDesc, $length1, $length2, $angle);

format STDOUT_TOP =
    Bond        Length A-B  Length B-C    Angle (B-A-C)
.

format STDOUT =
@|||||||||||||  @###.##     @###.##     @###.##
$bondDesc, $length1, $length2, $angle
.

my $atoms = $Documents{"urea.xsd"}->AsymmetricUnit->Atoms;

foreach my $atom (@$atoms) {
  #Get the bonded atoms
  my $atomBonds = $atom->Bonds;

  #Get number of distinct bonds
  my $bondCount = $atomBonds->Count;

  if($bondCount == 1) {
    next; #Ignore
  }
  elsif($bondCount > 1) {
    for(my $i = 0; $i < $bondCount; ++$i) {
      for(my $j = $i + 1; $j < $bondCount; ++$j) {
        #Atoms can be bonded in any order
        my $atom1 = $atomBonds->Item($i)->Atom1;
        if($atom->Name eq $atom1->Name) {
          $atom1 = $atomBonds->Item($i)->Atom2;
        }

        my $atom2 = $atomBonds->Item($j)->Atom1;
        if($atom->Name eq $atom2->Name) {
          $atom2 = $atomBonds->Item($j)->Atom2;
        }

        if(!($atom1->IsSymmetryParent && $atom2->IsSymmetryParent) &&
           !($atom1->Name eq $atom2->Name)) {
           #Ignore atoms that are not part of the symmetry parent unless,
           #they are images of the same atom
           next;
        }

        $bondDesc = $atom1->Name . "-" . $atom->Name . "-" . $atom2->Name;
        $length1  = $atomBonds->Item($i)->Length;
        $length2  = $atomBonds->Item($j)->Length;
        $angle    = rad2deg(CalcAngle($atom->XYZ, $atom1->XYZ, $atom2->XYZ));

        write; #Invoke Formatted Output
      }
    }
  }
}
