#!perl
#
# Purpose: Introduces substitutional disorder by replacing a percentage of atoms with another element. 
#          If this is a zeolite, Lowenstein's rule can be used to avoid too many substitutions.
#          Note. When calculating the number of atoms to change, this number is rounded up. This can
#          be changed by modifing the line containing $numToChange.

use strict;
use warnings;
use List::Util qw(shuffle); 	# Allows randomizing of the atoms in an array
use POSIX; 			# Allows use of ceil to round-up numbers
use MaterialsScript qw(:all);

# Lowenstein's rule applies to substitution of Aluminium in zeolite and states that no two
# Al atoms can share a common oxygen.

sub CheckLowenstein {
	my ($atom) = @_;
	
	#Rule states that Al-O-Al linkages are forbidden
	my $attached = $atom->AttachedAtoms;
	
	foreach my $firstAtom (@$attached) {
		#Only care about oxygens
		if($firstAtom->ElementSymbol eq "O") {
			my $oxyAttached = $firstAtom->AttachedAtoms;
			foreach my $secondAtom (@$oxyAttached) {
				if($secondAtom->ElementSymbol eq "Al") {
					
					return 0; #Lowenstein broken
					
				}
			}
		}
	}
	return 1; #Lowenstein satisfied
}

# Main routine to substitute atoms

sub SubstituteAtoms {

	my ($doc, $pcChange, $original, $new, $lowenstein)= @_;
	
	# Change the symmetry to P1
	
	$doc->MakeP1;
	
	srand;
	
	# Grab the collection of atoms in the unit cell
	
	my $atoms = $doc->UnitCell->Atoms;
	
	# Store the original atoms in an array
	
	my @originalAtoms;
	
	foreach my $atom (@$atoms) {
		if($atom->ElementSymbol eq "$original") {push(@originalAtoms, $atom);}
	}
	

	# Check to see if there are any atoms to modify
	
	if (scalar @originalAtoms == 0) {die "There are no atoms of element $original to modify\n";}

	# Randomize the atom positions in the array
	
	my @shuffledOriginalAtoms = shuffle(@originalAtoms);
	
	
	# calculate the number of atoms to change. The ceil commanound rounds up. If you wish to
	# round down, use int instead of ceil.
	
	my $numToChange = ceil (($pcChange/100) * scalar @shuffledOriginalAtoms);
	
	# Iterate through a counter and search for an atom to change
	
	for(my $i=0; $i<$numToChange; ++$i) {
		
		# Check to see if there are any atoms to modify
	
		if (scalar @shuffledOriginalAtoms == 0) {die "Modified $i atoms. There are no atoms of element $original to modify\n";}

		my $attempts = 0;
	
		# Look for an atom to modify
		
		my $atom=undef;
		
		while(!$atom && scalar @shuffledOriginalAtoms > 0) {
	
			my $index = int rand( scalar @shuffledOriginalAtoms);
	
			$atom = $shuffledOriginalAtoms[$index];
	
			# Check we haven't already changed it and it satisfies Lowenstein rule (optionally)
	
			if($atom->ElementSymbol eq "$new") {
	
				splice @shuffledOriginalAtoms, $index, 1; $atom = undef;
								
			} elsif ($lowenstein eq "Yes") { 
			
				if (!CheckLowenstein($atom)) {
	
					splice @shuffledOriginalAtoms, $index, 1; $atom = undef;	
				}	
			}
		
			++$attempts;
			
		}
	
		# Stop if atom is still undefined as this means Lowensteins rule has been broken
		
		if(!$atom) {
			print "Unable to satisfy Lowenstein's rule\n";
			
			# Count the number of atoms changed and report this.
			
			my $modAtoms = 0;
			
			foreach my $at (@{$doc->UnitCell->Atoms}) {
			
				if ($at->ElementSymbol eq "$new") { ++$modAtoms;}
			
			}
			
			print "$modAtoms atoms have been modified to $new element\n";	
			
			return;
		} else {
	
			#Change the element and display style
			$atom->ElementSymbol = "$new";
			$atom->Style         = "Ball and Stick";
	
		}
	
	}

}

# Specify the input document and percentage number of silicons to change

my $xsd = $Documents{"TON.xsd"};# Structure name
my $percentChange = 15;		# Percentage of atoms of original element to change
my $originalElement = "Si";	# Original element to change to new element
my $newElement = "Al";		# New element
my $obeyLowenstein = "Yes";	# Whether to obey Lowenstein's for zeolites

SubstituteAtoms($xsd, $percentChange,$originalElement, $newElement, $obeyLowenstein);
