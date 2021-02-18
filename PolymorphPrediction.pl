#!perl

use strict;
use MaterialsScript qw(:all);

sub DoGeometryOptimization($);

##################################################################
########       User defined parameters              ##############
##################################################################

# the name of the file containing the asymmetric unit
my $asymUnitFile = "ohd.xsd";

# Please see the Materials Studio Scripting help for more settings

# Energy server related settings
my @energyServerSettings = (
  CurrentForcefield     => 'Dreiding',
  ChargeAssignment      => 'Use current',
# KeepMotionGroupsRigid => 'true' # if set to true, rigid body minimization is done during geometry optimization
        );


# Settings specific to Polymorph
my @pmpSettings = (

# COARSE setting for testing purposes only, change to FINE for production runs
  Quality                    => "Coarse",
# MaxSteps                   => 50,  # maximum number of packing steps
# MaxClusters                => 2000,   # maximum number of clusters (not for Ultra-fine quality, where all clusters are used

# Reproducibility for testing purposes, comment out for production runs
  RandomNumberSeed           => 178234,

  @energyServerSettings
         );

# Should there be a clustering step after packing?
my $doPreclustering = 1;   # Yes
# my $doPreclustering = undef; # No

# Should the input asymmetric unit be optimized before starting a prediction?
# Only recommended if there is a single molecule in the document.
my $optimizeAsymUnit = 1;   # Yes
# my $optimizeAsymUnit = undef; # No

# List of space groups to be used.
# my @spaceGroups = (14 , 2, 19, 15, 4, 61, 33, 9, 60, 5 ); # 10 most frequent, using inernational table numbers
# my @spaceGroups = ("P21/C" , "P-1", "P212121", "C2/C", "P21"); # 5 most frequent, using space group symbols
my @spaceGroups = ("P212121" , "PBCA"); # the two space groups used in the tutorial

##################################################################
##################################################################

Modules->Polymorph->ChangeSettings(\@pmpSettings);

my $asymUnitDoc = $Documents{$asymUnitFile} or die "The document $asymUnitFile could not be loaded or found.";


my $seedName = $asymUnitDoc->Name;

# all results are added into this single study table.
my $resultTable = Documents->New("$seedName.std") or die "Unable to create a result study table $seedName.std";


print "\nPolymorph Prediction for $seedName \n\n";
print "Results are written into study table: $seedName.std \n\n";


# optimze the input asymmetric unit. Only use if the asymmetric unit contains a single fragment.
# Highly recommended, in particular if rigid body optimization is used.
if ( $optimizeAsymUnit )
{
        Modules->Forcite->ChangeSettings(\@energyServerSettings); #same energy server settings as Polymorph
        Modules->Forcite->GeometryOptimization->Run( $asymUnitDoc );

        print "\nInput asymmetric unit \"$seedName\" has been optimized.\n";
        printf "Total energy: : %7.3f \n" ,$asymUnitDoc->PotentialEnergy;
}


# the Predict cycle

foreach my $spaceGroup (@spaceGroups)
{

  print "\n ====Prediction for space group ID $spaceGroup ==== \n\n";

  my $packingResults = Modules->Polymorph->Packing->Run($asymUnitDoc, $spaceGroup);

  my $clusteredTrj;

  if ( $doPreclustering )
  {
        my $clusteringResults = Modules->Polymorph->Clustering->Run($packingResults->Trajectory);
        $clusteredTrj = $clusteringResults->Trajectory;
  }
  else
  {
        $clusteredTrj = $packingResults->Trajectory;
  }

  my $minTrj = DoGeometryOptimization( $clusteredTrj );

  my $clminResults = Modules->Polymorph->Clustering->Run($minTrj);

  Modules->Polymorph->Analysis->AddToStudyTable($clminResults->Trajectory , $resultTable);

  print "\n   Final trajectory document     : " . $clminResults->Trajectory->Name . "\n";
  print "   Number of predicted structures: " . $clminResults->Trajectory->NumFrames . "\n";
}

print "\n\n Polymorph script finished.\n\n";


########################################################################
#
#
# Performs a geometry optimization.
#
sub DoGeometryOptimization($) {

  my ($trj)=@_;

  my $results = Modules->Polymorph->GeometryOptimization->Run($trj);
  return $results->Trajectory;
}



