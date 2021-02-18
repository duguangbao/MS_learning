#!perl

use strict;
use MaterialsScript qw(:all);


sub AddRefineStructureDOF($);

#Load the settings file.
#The Reflex Refinement dialog has been used to configure the basic Refinement
#options to be applied to each phase. The settings have been saved using the 
#Save Settings ... options and entering "QPA" for the name of the settings.
Modules->Reflex->LoadSettings("QPA");


print   "\n\n================ QPA using the Rietveld method =======================\n";

my $phaseTable2 = Documents->New("StructureQPA.std");


#Add the first structure using the Refinement settings as saved in the QPA settings file.

my $struct1 =  $Documents{"Acetohexamide_FormA.xsd"};

AddRefineStructureDOF($struct1);

#turning off preferred refinement optimization in this particular case.
my $results = Modules->Reflex->PowderQPA->AddPhase($phaseTable2, $struct1 ,Settings( PreferredOrientation => "None"));
my $phaseRowS1 = $results->RowIndex;


#Phase B

my $struct2 =  $Documents{"Acetohexamide_FormB.xsd"};

AddRefineStructureDOF($struct2);                            
                              
$results = Modules->Reflex->PowderQPA->AddPhase($phaseTable2, $struct2);
my $phaseRowS2 = $results->RowIndex;


#corrundum

my $struct3 =  $Documents{"corundum.xsd"};

#Delete all bonds as those prevent automatic definition of structural freedom.
$struct3->AsymmetricUnit->Bonds->Delete();
AddRefineStructureDOF($struct3);

#corundum has atoms on special positions - during lattice change, fractional coordinates need to be fixed                                                                                      
$results = Modules->Reflex->PowderQPA->AddPhase($phaseTable2, $struct3, Settings( FixFractionalCoordinates => "Yes"));
my $phaseRowS3 = $results->RowIndex;


#==== Modify the settings for the QPA calculation if necessary

Modules->Reflex->ChangeSettings(
                     Settings(Convergence               => "Coarse", 
                              AddPhaseDiffractionCharts => "No",
                              TwoThetaMin               =>  5,
                              TwoThetaMax               =>  80
                              ));



#==== Load the experimental data and run the QPA calculation. If multiple mixture pattern with the same phases
#==== are to be decomposed, just start a loop here.

# my $expDoc = Documents=>Import("C:\MyPatterns\Acetohexamide_Mixture.3cam"); #use the import method if you need to import
                                                                              #data from outside the current project
my $expDoc = $Documents{"Acetohexamide_Mixture.xcd"};


my $QPAResults = Modules->Reflex->PowderQPA->Run($phaseTable2, $expDoc); 


printf   "Rwp for decomposition of %s: %5.1f \n", $expDoc->Name , $QPAResults->Rwp ;

my $weightA = $phaseTable2->Cell( $phaseRowS1, 'Weight %');
my $weightB = $phaseTable2->Cell( $phaseRowS2, 'Weight %');
my $weightC = $phaseTable2->Cell( $phaseRowS3, 'Weight %');

printf   "\nContent of Form A  : %6.2f \n", $weightA ;
printf     "Content of Form B  : %6.2f \n", $weightB ;
printf     "Content of corundum: %6.2f \n", $weightC ;

print "\nQPA calculation finished.\n\n";


#======================================
#
# Adds structural DOF and sets refinement flags, no torsion refinement

sub AddRefineStructureDOF($)
{

	my ($struct) = @_;

	Modules->Reflex->Preparation->SetRefineLattice( $struct, "Yes" );
	                              
	Modules->Reflex->Preparation->AssignStructuralDegreesOfFreedom(  $struct );
	Modules->Reflex->Preparation->SetRefineTorsions(  $struct, "No" );
	Modules->Reflex->Preparation->SetRefineMotionGroups(  $struct, "Yes" ) ;

}
