#!perl
#
# Purpose: Demonstrate processing each atomistic document in a project in turn

use strict;
use warnings;
use MaterialsScript qw(:all);

sub CountAtoms {
    my ($doc) = @_;
    
    my $docname = $doc->Name.".xsd" ;
    
    my %elementStats;
    
    foreach my $atom (@{$doc->UnitCell->Atoms}) {
        $elementStats{$atom->ElementSymbol} += 1;    
    }
    
    print "Atoms stats for $docname:\n";
    
    foreach my $element (sort keys %elementStats) {
    	print "\t $element $elementStats{$element}\n";
    }
    
    $doc->Close;
}

#Filter out only the .xsd files in the Documents collection
foreach my $key (keys %Documents) {
    my $doc = $Documents{$key}; 
    
    if ($doc->Type eq "3DAtomistic" ) { 
       CountAtoms($doc);
    }
}
