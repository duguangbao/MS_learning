#!perl
#
# Purpose: Replace 10% of the silicon atoms in a zeolite with aluminium.
#          A sodium atom is then bonded to any aluminium atoms introduced
#          into the structure.

use strict;
use warnings;
use MaterialsScript qw(:all);

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

sub BuildZeolite {
    my ($docname) = @_;
    
    my $doc = $Documents{$docname};
    $doc->MakeP1;
    
    srand;
    
    my $atoms = $doc->UnitCell->Atoms;
    
    my @siAtoms;
    
    #First find all the Silicon atoms
    foreach my $atom (@$atoms) {
        if($atom->ElementSymbol eq "Si") {
            push(@siAtoms, $atom);
        }
    }
    
    #Modify 10% of the Silicon Atoms
    my $numToChange = int 0.1 * @siAtoms;
    print "\tThere are ", scalar @siAtoms, " silicon atoms. Modifying $numToChange.\n";
    for(my $i=0; $i<$numToChange; ++$i) {
        my $atom;
    
        #Look for an atom to modify
        while(!$atom && scalar @siAtoms > 0) {
            my $index = int rand(@siAtoms);
            #print "Checking Atom at $index\n";
    
            $atom = $siAtoms[$index];
    
            #Check we haven't already changed it and it satisfies Lowenstein rule
            if($atom->ElementSymbol eq "Al" || !CheckLowenstein($atom) ) {
                splice @siAtoms, $index, 1;
                
                $atom = undef;
            }
        }
        
        if(!$atom) {
            print "\tUnable to Satisfy Lowenstein's rule in input structure";
            
            $doc->Close;
            
            return;
        }
    
        #Change the element and display style
        $atom->ElementSymbol = "Al";
        $atom->Style         = "Ball and Stick";
    
        my $pos = $atom->XYZ;
    
        #Create a new Sodium somewhere near Al
        my $sodium = $atoms->CreateAtom("Na", Point( X=> $pos->X + rand(2),
                                                     Y=> $pos->Y + rand(2),
                                                     Z=> $pos->Z + rand(2) ));
      
        $sodium->Style = "Ball and Stick";
        $atoms->CreateBond($atom, $sodium, "Single");
    }
    
    Modules->Forcite->GeometryOptimization->Run($doc, Settings(MaxIterations => 1000));    
}

BuildZeolite("TON.xsd");
