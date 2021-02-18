#!perl
#####################################################################################
# Description: Generates coarse-grained interactions for a bead trajectory
#
# Authors: Neil Spenley, Reinier Akkermans (Biovia)
# Version: 3.2 (July 2018)
# Modules: Materials Visualizer, Mesocite
#
# This script takes as input a bead trajectory (xtd), typically a coarse-grained
# atomistic trajectory obtained from a Molecular Dynamics simulation at constant 
# temperature using Forcite.
#
# The script generates a coarse-grained force field which reproduces the statistics
# of the input trajectory. This is done using an iterative procedure in which dynamics
# is run repeatedly and the result is used to update the force field.
#
# For each bond/angle/nonbond between distinct types an interaction is introduced, based
# on a tabulated potential. The script calculates the distribution function for each interaction
# in the trajectory, converts this to a potential of mean force, and adds it to the forcefield as a
# tabulated potential. Then dynamics is run with the new forcefield, the result is compared
# with the reference trajectory, and the potential is modified. 
#
# Required input:
# 1. Coarse-grained trajectory (coarse-grained atomistic trajectory)
# 2. Temperature
#
# Output documents:
# 1. A sequence of force field documents containing the force field parameters
# 2. Study table containing all fitting data 
#####################################################################################
use strict;
use warnings;
use MaterialsScript qw(:all);
 
{
################################################################################
#  BEGIN USER INPUT                                                            #
################################################################################
# the trajectory to be analysed
my $referenceName = "Input"; # name of input reference trajectory
my $trialName = "Fitting";       # root name for documents created by this script


# script settings
my $settings = Settings
(
	UseBondStretch => "Yes",
	LengthDistributionBinWidth => 0.01, # binwidth for length distribution (in A)

	UseAngleBend => "Yes",
	AngleDistributionBinWidth => 0.5, #  binwidth for angle distribution (in deg)

	UseTorsion => "No",
	TorsionDistributionBinWidth => 5, #  binwidth for torsion distribution (in deg)
	
	UseVanDerWaals => "Yes",
	RDFBinWidth => 0.1, #  bin width for radial distribution distribution (in A)
	RDFCutoff => 15, # maximum distance for the radial distribution distribution (in A)
	VanDerWaalsCombinationRule => "None", # if "None" cross-terms are also evaluated
	
	Temperature => 298, # temperature at which the trajectory was created (in K)
	
	DistributionCutoff => 0.005, # ignore points less than this fraction of maximum probability
	
	NumberOfIterations => 10,
	NumberOfSteps => 2000000,    # number of timesteps in each dynamics run
	
	RelaxationFactor => 1,       # controls the update of the potentials in each iteration:
	                             # new potential = old potential + (relaxation factor)*(reference pmf - last trial pmf)
	NumZeroesUntilCutoff => 4,  # candidate nonbond potentials can be truncated after the nth zero to avoid long range oscillations
	                             # For short range repulsion and one attractive well, set to 2
	
	StartIteration => 1,       # Normally 1. Can use to restart a run from a later iteration
	PartialForcefield => "",   # Optional existing forcefield with some interactions. New interactions will be added to a copy of this
);
################################################################################
#  END USER  INPUT                                                             #
################################################################################

# prepare the structure for analysis
my $referencedoc = $Documents{"$referenceName.xtd"};
my @types = PrepareStructure( $referencedoc );

# create a hash of the settings array for easy access
my %s = @$settings;
	
my $startIteration = $s{"StartIteration"} ;
if ( $startIteration == 1 )
{
	my $doc = Documents->New("$trialName" . "0.xtd");
	$doc->CopyFrom( $referencedoc );
	RunAnalysis($trialName, 0, $doc, undef, $settings, @types);
}

for my $iteration ($startIteration..$s{"NumberOfIterations"})
{
	my $temperature = $s{"Temperature"};
	
	my $ff = CreateForcefield($trialName, $referencedoc, $settings, $iteration, @types);
	my $doc = Documents->New("$trialName$iteration.xsd");
	$doc->CopyFrom( $referencedoc );
	PrepareStructure($doc);
	my $results = Modules->Mesocite->Dynamics->Run($doc, Settings(
		CurrentForcefield => '/' . $ff->Name, 
		Temperature => $temperature,
		Ensemble3D => "NVT", 
		NumberOfSteps => $s{"NumberOfSteps"}, 
		TrajectoryFrequency => 1000, 
		TimeStep => 1,
		Thermostat => "NHL"));

	if ( abs( $results->Temperature - $temperature ) > 5 ) { die "Temperature out of control\n"; }
	my $trajectorydoc = $Documents{"$trialName$iteration.xtd"};
	RunAnalysis($trialName, $iteration, $trajectorydoc, $ff, $settings);
}
printf "Time taken: %d seconds", time()-$^T;
}

################################################################################
#  Main routine
################################################################################
sub RunAnalysis
{
	my ($trialName, $iteration, $doc,$trialFf,$settings, @types) = @_;

	# create a hash of the settings array for easy access
	my %s = @$settings;
	
	# some constants
	my $kT = $s{Temperature}*8.314/4.184/1000; # kcal/mol
	my $analysis = Modules->Mesocite->Analysis;
	
	# prepare the structure for analysis
	PrepareStructure($doc);

	my $std = FindOrCreateStudyTable( $trialName );
	
	# fit the bonds, if used
	if ($s{"UseBondStretch"} eq "Yes")
	{
		$analysis->ChangeSettings([LengthDistributionBinWidth => $s{LengthDistributionBinWidth}]);
		
		my @typeSequences = PrepareBonds($doc);

		foreach my $typeSequence (@typeSequences)
		{			
			# analyze the trajectory and add the distribution to the table
			LengthDistribution($std, $iteration, $doc, $trialFf, $typeSequence, $kT, \%s );
		}
	}

	# fit the angles, if used
	if ($s{"UseAngleBend"} eq "Yes")
	{
		$analysis->ChangeSettings([AngleDistributionBinWidth => $s{AngleDistributionBinWidth}]);
	
		my @typeSequences = PrepareAngles($doc);	

		foreach my $typeSequence (@typeSequences) 
		{			
			# analyze the trajectory and add the distribution to the table
			AngleDistribution($std, $iteration, $doc, $trialFf, $typeSequence, $kT, \%s);
		}
	}
	
	# fit the torsions, if used
	if ($s{"UseTorsion"} eq "Yes")
	{
		$analysis->ChangeSettings([TorsionDistributionBinWidth => $s{TorsionDistributionBinWidth}]);
		
		my @typeSequences = PrepareTorsions($doc);

		foreach my $typeSequence (@typeSequences)
		{			
			# analyze the trajectory and add the distribution to the table
			TorsionDistribution($std, $iteration, $doc, $trialFf, $typeSequence, $kT, \%s );
		}
	}
	
	# fit the non-bonds, if used
	if ($s{"UseVanDerWaals"} eq "Yes")
	{
		$analysis->ChangeSettings([
			RDFBinWidth => $s{RDFBinWidth},	RDFCutoff => $s{RDFCutoff},
			RDFComputeMolecularComponents => "Yes"]);
		
		my $rule = $s{VanDerWaalsCombinationRule};
		my $useCrossterms = $rule eq "None";

		my @typeSequences = PrepareNonBonds($doc,$useCrossterms);

		foreach my $typeSequence (@typeSequences) 
		{	
			# analyze the trajectory and add the distribution to the table
			NonBondDistribution($std, $iteration, $doc, $trialFf, $typeSequence, $kT, \%s);
		}
	}
	
	$std->Save;
}
#########################################################################
sub CreateForcefield
{
	my ($trialName, $doc, $settings, $iteration, @types) = @_;

	# create a hash of the settings array for easy access
	my %s = @$settings;
	
	# some constants
	my $kT = $s{Temperature}*8.314/4.184/1000; # kcal/mol
	my $analysis = Modules->Mesocite->Analysis;

	# create a new force field
	my $off = Documents->New($trialName . $iteration . ".off");

	if (defined $s{"PartialForcefield"} && $s{"PartialForcefield"} ne "" )
	{
		$off->CopyFrom( $Documents{ $s{"PartialForcefield"} } );
	}
	else
	{
		# set up the types in the force field
		PrepareForcefield($off, \@types, $s{UseBondStretch},$s{UseAngleBend},$s{UseTorsion},$s{UseVanDerWaals});
	}
	
	# the study table output document
	my $std = $Documents{$trialName . ".std"};
	
	PrepareStructure($doc);
	
	# fit the bonds, if used
	if($s{"UseBondStretch"} eq "Yes")
	{	
		my @typeSequences = PrepareBonds($doc);

		foreach my $typeSequence (@typeSequences)
		{			
			my $name = "bond " .  $typeSequence;
			my $potential = EstimatePotential($std, $name, $iteration-1, $s{RelaxationFactor}, $s{DistributionCutoff} );
			FillGapsInPotential( $potential, "QuadraticAtEnds" );
			$off->CreateTerm("Bond Stretch", $typeSequence, "Tabulated", [ R => @{$potential->[0]}, E => @{$potential->[1]} ]);
		}
	}

	# fit the angles, if used
	if($s{"UseAngleBend"} eq "Yes")
	{
		my @typeSequences = PrepareAngles($doc);	

		foreach my $typeSequence (@typeSequences) 
		{			
			# fit the distribution data
			my $name = "angle " . $typeSequence;
			my $potential = EstimatePotential($std, $name, $iteration-1, $s{RelaxationFactor}, $s{DistributionCutoff} );
			FillGapsInPotential( $potential, "QuadraticAtEnds" );
			$off->CreateTerm("Angle Bend", $typeSequence, "Tabulated", [ T => @{$potential->[0]}, E => @{$potential->[1]} ]);
		}
	}
	
	# fit the torsions, if used
	if($s{"UseTorsion"} eq "Yes")
	{
		my @typeSequences = PrepareTorsions($doc);	

		foreach my $typeSequence (@typeSequences) 
		{			
			# fit the distribution data
			my $name = "torsion " . $typeSequence;
			my $potential = EstimatePotential($std, $name, $iteration-1, $s{RelaxationFactor}, $s{DistributionCutoff} );
			FillGapsInPotential( $potential, "Periodic" );
			$off->CreateTerm("Torsion", $typeSequence, "Tabulated", [ Phi => @{$potential->[0]}, E => @{$potential->[1]} ]);
		}
	}
	
	# fit the non-bonds, if used
	if($s{"UseVanDerWaals"} eq "Yes")
	{
		my $rule = $s{VanDerWaalsCombinationRule};
		my $useCrossterms = $rule eq "None";
		$off->VanDerWaalsCombinationRule = $rule if !$useCrossterms;

		my @typeSequences = PrepareNonBonds($doc,$useCrossterms);

		foreach my $typeSequence (@typeSequences) 
		{	
			my $name = "non-bond " . $typeSequence;
			my $potential = EstimatePotential($std, $name, $iteration-1, $s{RelaxationFactor}, $s{DistributionCutoff} );
			FillGapsInPotential( $potential, "TruncateAtLargeX" );
			my $truncatedPotential = ZeroPotentialAfterNthZero( $potential, $s{"NumZeroesUntilCutoff"} );
			my $rData = $truncatedPotential->[0];
			my $EData = $truncatedPotential->[1];
			$off->CreateTerm("van der Waals", $typeSequence, "Tabulated", [ R => @$rData, E => @$EData ]);
		}
	}
	
	$std->Save;
	
	return $off;
}
#########################################################################
# Find all unique connectors and create sets
sub PrepareBonds
{
	my ($doc) = @_;
	
	my @names = ();
	my $connectors = $doc->UnitCell->BeadConnectors;
	foreach my $connector (@$connectors) 
	{
		my $bead1 = $connector->Bead1;
		my $bead2 = $connector->Bead2;
		if ($bead1->ForcefieldType eq "" || $bead2->ForcefieldType eq "" ) { die "Forcefield types not defined" };
		# create a distance monitor for this connector
		my $distance = $doc->CreateDistance([$bead1,$bead2]);
		
		# name the monitor by sequence; use a fixed order
		my @types;
		push @types,$bead1->ForcefieldType;
		push @types,$bead2->ForcefieldType;
		my @typesSorted = sort(@types);
		my $name = $typesSorted[0].",".$typesSorted[1];
		$distance->Name = $name;
		push @names, $name;
	}

	my @namesFiltered = CreateSetForEachName($doc->UnitCell->Distances, \@names, $doc);
	return @namesFiltered;	
}
#########################################################################
# Find all unique angles and create sets
sub PrepareAngles
{
	my ($doc) = @_;
	
	my @names = ();
	my $beads = $doc->UnitCell->Beads;
	foreach my $bead (@$beads) 
	{
		my $type2 = $bead->ForcefieldType;
	
		# loop over all pairs of neighbours
		my $attachedBeads = $bead->AttachedBeads;
		my $count = $attachedBeads->Count;
		for(my $i = 0; $i < $count-1; $i++)
		{
			my $bead1 = $attachedBeads->Item($i);
			my $type1 = $bead1->ForcefieldType;
			
			for(my $j = $i+1; $j < $count; $j++)
			{
				my $bead3 = $attachedBeads->Item($j);
				my $type3 = $bead3->ForcefieldType;

				# create an angle monitor for this angle
				my $angle = $doc->CreateAngle([$bead1,$bead,$bead3]);
				
				# name the monitor by sequence; use a fixed order
				my @types;
				push @types,$type1;
				push @types,$type3;
				my @typesSorted = sort(@types);
				my $name = $typesSorted[0].",".$type2.",".$typesSorted[1];
				$angle->Name = $name;
				push @names, $name
			}
		}
	}

	my @namesFiltered = CreateSetForEachName($doc->UnitCell->Angles, \@names, $doc);
	return @namesFiltered;
}

#########################################################################
# Find all unique torsions and create sets
sub PrepareTorsions
{
	my ($doc) = @_;
	
	AssignBeadNames($doc);
	
	my @names = ();
	my $connectors = $doc->UnitCell->BeadConnectors;
	foreach my $connector (@$connectors) 
	{
		my $bead1 = $connector->Bead1;
		my $bead2 = $connector->Bead2;
		if ($bead1->ForcefieldType eq "" || $bead2->ForcefieldType eq "" ) { die "Forcefield types not defined" };
		my @namesForThisConnector = CreateAllTorsionsForBeadPair($bead1, $bead2, $doc);
		push @names, @namesForThisConnector;
	}

	my @namesFiltered = CreateSetForEachName($doc->UnitCell->Torsions, \@names, $doc);
	return @namesFiltered;
}
#########################################################################
sub CreateAllTorsionsForBeadPair
{
	my ($beadJ, $beadK, $doc) = @_;
	
	my $attachedBeadsToJ = $beadJ->AttachedBeads;
	my $attachedBeadsToK = $beadK->AttachedBeads;
	
	my @names;
	for my $beadI (@$attachedBeadsToJ)
	{	
		for my $beadL (@$attachedBeadsToK)
		{
		    # Avoid the case where bead I is the same object as bead K (or L is the same as J)
			if (($beadI->Name ne $beadK->Name) && ($beadL->Name ne $beadJ->Name))
			{
				my $torsion = $doc->CreateTorsion([$beadI,$beadJ,$beadK,$beadL]);
			
				# name the monitor by sequence; use a fixed order
				my @types = ( $beadI->ForcefieldType, $beadJ->ForcefieldType, $beadK->ForcefieldType, $beadL->ForcefieldType );
				if ( $types[1] gt $types[2] )
				{
					@types = reverse (@types);
				}

				my $name = join( ',', @types );
				$torsion->Name = $name;
				push @names, $name;
			}
		}
	}
	
	return @names;
}
#########################################################################
sub AssignBeadNames
{
	my ($doc) = @_;
	my $beadIndex = 0;
	for my $bead (@{$doc->UnitCell->Beads})
	{
		$bead->Name = "$beadIndex";
		$beadIndex++;
	}
}
#########################################################################
sub CreateSetForEachName
{
	my ($listOfObjects, $names, $doc) = @_;
	
	# filter for unique ones
	my %hash = ();
	my @uniqueNames = grep(!$hash{$_}++, @$names);

	foreach my $name (@uniqueNames) 
	{	 
		my @objectsWithSameName = ();
		foreach my $object (@$listOfObjects)
		{
			if($object->Name eq $name)
			{
				push @objectsWithSameName, $object;
			}
		}
		$doc->CreateSet($name, \@objectsWithSameName); 
	}

	return @uniqueNames;	
}
#########################################################################
# Create sets for each bead type
sub PrepareNonBonds
{
	my ($doc,$useCrossTerms) = @_;
	
	# find all force field types
	my @types = ();
	my $beads = $doc->UnitCell->Beads;
	foreach my $bead (@$beads) 
	{
		my $type = $bead->ForcefieldType;
		push @types, $type;
	}
	# filter for unique ones
	my %hash = ();
	my @typesFiltered = grep(!$hash{$_}++, @types);
	
	# create sets for each unique name
	foreach my $type (@typesFiltered) 
	{
		# filter all distance for the one we want
		my $set = $doc->CreateSet($type,[]);
		
		foreach my $bead (@$beads)
		{
			$set->Add($bead) if $bead->ForcefieldType eq $type
		}
	}

	# the self terms
	my @typeSequences = ();
	my $count = scalar(@typesFiltered);	
	for(my $i = 0; $i < $count; $i++)
	{
		my $type = $typesFiltered[$i];
		push @typeSequences, $type . "," . $type;
	}

	# the cross terms, if used
	if($useCrossTerms)
	{
		for(my $i = 0; $i < $count; $i++)
		{
			my $type1 = $typesFiltered[$i];
			for(my $j = $i+1; $j < $count; $j++)
			{
				push @typeSequences, $type1 . "," . $typesFiltered[$j];
			}
		}
	}
		
	return @typeSequences;	
}

################################################################################
# create new columns, add labels and x-data
sub FindOrCreateStudyTable
{
	my ($name) = @_;
	
	my $std;
	eval
	{
		$std = $Documents{$name . ".std"};
	};
	
	if (!defined $std)
	{
		$std = Documents->New($name . ".std");
		TrimSheet($std,0,5);
		$std->Title = "Summary";
		my $col = 0;
		$std->ColumnHeading($col++) = "Structure";
		$std->ColumnHeading($col++) = "Name";
		$std->ColumnHeading($col++) = "Distribution";
	}
	
	return $std;
}

####################################################################################
# Modifies the number of rows and columns of a sheet
sub TrimSheet
{
	my ($sheet,$nRow,$nCol) = @_;

	my $nRowDoc = $sheet->RowCount;
	if($nRowDoc > $nRow)
	{
		for(my $iRow = $nRowDoc-1; $iRow >= $nRow; $iRow--)
		{
			$sheet->DeleteRow($iRow);
		}
	}
	else
	{
		for (my $iRow = $nRowDoc; $iRow < $nRow; $iRow++)
		{
			$sheet->InsertRow($iRow);
		}		
	}

	my $nColDoc = $sheet->ColumnCount;
	if($nColDoc > $nCol)
	{
		for(my $iCol = $nColDoc-1; $iCol >= $nCol; $iCol--)
		{
			$sheet->DeleteColumn($iCol);
		}
	}
	else
	{
		for (my $iCol = $nColDoc; $iCol < $nCol; $iCol++)
		{
			$sheet->InsertColumn($iCol);
		}		
	}
}
#############################################################################
# set up the interactions, types and description
sub PrepareForcefield
{
	my ($off,$types,$useBondStretch,$useAngleBend,$useTorsion,$useVanDerWaals) = @_;

	$off->Description .= "\nCoarse-grained force field for " . $off->Name ."\n\n";

	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
	my $timeStamp = sprintf ("%4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);
	$off->Description .= "Created: $timeStamp\n\n";
	
	# create the types
	foreach my $type (@$types) 
	{
		$off->CreateType($type);
	}
	
	$off->UseBondStretch = $useBondStretch;
	$off->UseAngleBend = $useAngleBend;
	$off->UseTorsion = $useTorsion;
	$off->UseVanDerWaals = $useVanDerWaals;	
}
################################################################################
# ensures that each bead has a forcefield type assigned, creating new ones
# if undefined; returns a list of unique forcefield types
sub PrepareStructure
{
	my ($doc) = @_;
	
	if (!defined $doc) { return; }
	
	# collect names of all beads
	my $beads = $doc->UnitCell->Beads;
	my @names = ();
	foreach my $bead (@$beads) 
	{
		# if no force field types are defined, derive from bead type name
		my $name = $bead->ForcefieldType;
		if(length($name) <= 0)
		{
			my $beadTypeName = $bead->BeadTypeName;
					
			# the force field document restricts names to 5 characters
			$name = substr $beadTypeName, 0, 4;
			$bead->ForcefieldType = $name;
		}
	
		push @names, $name;
	}

	# filter for unique ones
	my %hash =();
	my @namesFiltered = grep(!$hash{$_}++, @names);
	
	return @namesFiltered;
}

################################################################################
# Calculate the distribution function for a given set
sub LengthDistribution
{
	my ($std,$iteration, $trialDoc,$ff,$type,$kT,$s) = @_;
	
	my $name = "bond " . $type;
	my $sheet = GetSheet($std, $name);
	if (!defined $sheet)
	{
		$sheet = AddNewSheetToStudyTable  ($std, $name, "r (A)", 0, $s->{"LengthDistributionBinWidth"}, 20 / $s->{"LengthDistributionBinWidth"} ); 
	}
	
	my ( $potentialColumnHeading, $pColumnHeading, $pmfColumnHeading ) = ColumnHeadingsForIteration( $iteration );

	CopyDataToResultsTable($sheet, $potentialColumnHeading, GetPotentialFromForcefield($ff, "Bond Stretch", $type, "R") ) if defined $ff;	
	
	my ($trialStdA, $trialChart) = GetDistributionAnalysis($trialDoc, "LengthDistribution", "LengthDistributionChart", Settings(LengthDistributionSetA => $type)); 
	if (!defined $trialStdA) { die "Could not calculate distribution for $name \n"; };
	AddDistributionAndPmfToResultsTable ($sheet,  $trialStdA, $kT, $pColumnHeading, $pmfColumnHeading,  "Distance", "P", sub{ my ($r) = @_; return ($r*$r); } );
	AddAnalysisResultToSummarySheet( $std, $type, $name, $trialDoc, $trialChart );
	CleanUpUnusedDocuments($trialStdA, $trialChart);
}
################################################################################
# Calculate the distribution function for a given set
sub AngleDistribution
{
 	use Math::Trig;
	my ($std,$iteration, $trialDoc,$ff,$type,$kT,$s) = @_;
	
	my $name = "angle " . $type;
	my $sheet = GetSheet($std, $name);
	if (!defined $sheet)
	{
		$sheet = AddNewSheetToStudyTable  ($std, $name, "angle (deg)", 0, $s->{"AngleDistributionBinWidth"}, 180 / $s->{"AngleDistributionBinWidth"} ); 
	}
	
	my ( $potentialColumnHeading, $pColumnHeading, $pmfColumnHeading ) = ColumnHeadingsForIteration( $iteration );
	
	CopyDataToResultsTable($sheet, $potentialColumnHeading, GetPotentialFromForcefield($ff, "Angle Bend", $type, "T" ) ) if defined $ff;

	my ($trialStdA, $trialChart) = GetDistributionAnalysis($trialDoc, "AngleDistribution", "AngleDistributionChart", Settings(AngleDistributionSetA => $type));
	if (!defined $trialStdA) { die "Could not calculate distribution for $name \n"; };
	AddDistributionAndPmfToResultsTable  ($sheet,  $trialStdA, $kT, $pColumnHeading, $pmfColumnHeading,  "Angle", "P", sub{ my ($angleDeg) = @_; return (sin($angleDeg*pi/180.)); } );
	AddAnalysisResultToSummarySheet( $std, $type, $name, $trialDoc, $trialChart );
	CleanUpUnusedDocuments($trialStdA, $trialChart);
}
################################################################################
# Calculate the distribution function for a given set
sub TorsionDistribution
{
 	use Math::Trig;
	my ($std,$iteration, $trialDoc,$ff,$type,$kT,$s) = @_;
	
	my $name = "torsion " . $type;
	my $sheet = GetSheet($std, $name);
	if (!defined $sheet)
	{
		$sheet = AddNewSheetToStudyTable  ($std, $name, "angle (deg)", -180, $s->{"TorsionDistributionBinWidth"}, 1 + 360 / $s->{"TorsionDistributionBinWidth"} ); 
	}
	
	my ( $potentialColumnHeading, $pColumnHeading, $pmfColumnHeading ) = ColumnHeadingsForIteration( $iteration );
	
	CopyDataToResultsTable($sheet, $potentialColumnHeading, GetPotentialFromForcefield($ff, "Torsion", $type, "Phi" ) ) if defined $ff;

	my ($trialStdA, $trialChart) = GetDistributionAnalysis($trialDoc, "TorsionDistribution", "TorsionDistributionChart", Settings(TorsionDistributionSetA => $type));
	if (!defined $trialStdA) { die "Could not calculate distribution for $name \n"; };
	AddDistributionAndPmfToResultsTable  ($sheet,  $trialStdA, $kT, $pColumnHeading, $pmfColumnHeading,  "Angle", "P", sub{ return 1; } );
	AddAnalysisResultToSummarySheet( $std, $type, $name, $trialDoc, $trialChart );
	CleanUpUnusedDocuments($trialStdA, $trialChart);
}
################################################################################
# Calculate the distribution function for a given set
sub NonBondDistribution
{
	my ($std,$iteration, $trialDoc,$ff,$type,$kT,$s) = @_;

	my $name = "non-bond " . $type;
	my $sheet = GetSheet($std, $name);
	if (!defined $sheet)
	{
		$sheet = AddNewSheetToStudyTable ($std, $name, "r (A)", 0.5*$s->{"RDFBinWidth"}, $s->{"RDFBinWidth"}, $s->{"RDFCutoff"} / $s->{"RDFBinWidth"} );
	}
	
	my ( $potentialColumnHeading, $pColumnHeading, $pmfColumnHeading ) = ColumnHeadingsForIteration( $iteration );
	
	my @types = split(',', $type);
	my $sheetName = "Inter";
	$sheetName .= " (" . $types[0] . " - " . $types[1] . ")" if $types[0] ne $types[1];

	CopyDataToResultsTable($sheet, $potentialColumnHeading, GetPotentialFromForcefield($ff, "van der Waals", $type, "R" ) ) if defined $ff;
		
	my ($trialStdA, $trialChart) = GetDistributionAnalysis($trialDoc, "RadialDistributionFunction", "RDFChart", Settings(RDFSetA => $types[0], RDFSetB => $types[1]) ); 
	if (!defined $trialStdA) { die "Could not calculate distribution for $name \n"; };
	my $trialSheet = GetSheet($trialStdA, $sheetName);
	AddDistributionAndPmfToResultsTable  ($sheet,  $trialSheet,     $kT, $pColumnHeading, $pmfColumnHeading,  "r", "g(r)", sub{ return 1; } );
	AddAnalysisResultToSummarySheet( $std, $type, $name, $trialDoc, $trialChart );
	CleanUpUnusedDocuments($trialStdA, $trialChart);
}
################################################################################
sub GetDistributionAnalysis
{
	my ($doc, $distributionType, $outputChartName,  $settings) = @_;
	if (defined $doc)
	{
		my $analysis = eval 'Modules->Mesocite->Analysis->' . $distributionType . '($doc,$settings)';
		return ( $analysis->{ "${outputChartName}AsStudyTable" }, $analysis->{ "${outputChartName}" } );
	}
	else
	{
		return (undef,undef);
	}
}
################################################################################
sub CleanUpUnusedDocuments
{
	for my $doc (@_)
	{
		if (defined $doc)
		{
			$doc->Delete;
		}
	}
}
################################################################################
sub AddNewSheetToStudyTable 
{
    my ($std, $name, $xColumnHeading, $xBegin, $deltaX, $rowCount ) = @_;
	my $sheet = $std->InsertSheet;
	$sheet->Title = $name;
	TrimSheet($sheet,$rowCount,1);
	$sheet->ColumnHeading(0) = $xColumnHeading;
	for(my $row = 0; $row < $rowCount; $row++)
	{
		my $x = $xBegin + $row*$deltaX;
		$sheet->Cell($row,0) = $x;
	}
	return $sheet;
}
################################################################################
sub CalculatePMF
{
    my ($kT, $distribution, $jacobianCorrectionFunction ) = @_;
	
	my $xData = $distribution->[0];
	my $probabilityData = $distribution->[1];
	my @pmfData;
	for my $row (0..$#$xData)
	{
		my $x = $xData->[$row];
		my $Prob = $probabilityData->[$row];
        my $jacobianCorrection = &$jacobianCorrectionFunction($x);
		if ( $Prob > 0 &&  $jacobianCorrection > 0)
		{
			push @pmfData, -$kT*log($Prob / $jacobianCorrection );
		}
		else
		{
			push @pmfData, undef;
		}
	}
	return [ [@$xData], \@pmfData ];
}
################################################################################
sub GetDistributionFromAnalysisStudyTable
{
    my ( $stdA, $xVariableInAnalysisOutput, $yVariableInAnalysisOutput ) = @_;
	
	my $rowCount = $stdA->RowCount;
	my @xdata;
	my @distributionData;
	for (my $row = 0; $row < $rowCount; $row++)
	{
		my $r = $stdA->Cell($row,$xVariableInAnalysisOutput);
		push @xdata, $r;
		my $Prob = $stdA->Cell($row,$yVariableInAnalysisOutput);
		push @distributionData, $Prob;
	}
	return [ \@xdata, \@distributionData ];
}
################################################################################
sub AddDistributionAndPmfToResultsTable
{
    my ($sheet, $stdA, $kT, $xColumnHeading, $yColumnHeading, $xVariableInAnalysisOutput, $yVariableInAnalysisOutput, $jacobianCorrectionFunction ) = @_;

	my $distribution = GetDistributionFromAnalysisStudyTable( $stdA, $xVariableInAnalysisOutput, $yVariableInAnalysisOutput );
	my $pmf = CalculatePMF( $kT, $distribution, $jacobianCorrectionFunction );
	CopyDataToResultsTable($sheet, $xColumnHeading, $distribution );
	CopyDataToResultsTable($sheet, $yColumnHeading, $pmf );
}
################################################################################
sub CopyDataToResultsTable
{
    my ($outputSheet, $columnHeading, $data) = @_;
	# $data is pair of arrays: [ [x1,x2,x3,...], [y1,y2,y3,...] ]
	# Output sheet has x values already in column 0

	my $dataColumn = $outputSheet->InsertColumn($outputSheet->ColumnCount, $columnHeading );	
	my $xData = $data->[0];
	my $yData = $data->[1];
	
	my $tol = 1e-8;
	
	my $inputRow = 0;
	my $outputRow = 0;
	while ( $inputRow < scalar(@$xData) && $outputRow < $outputSheet->RowCount )
	{
		my $inputx = $xData->[$inputRow];
		
		if ( abs( $inputx - $outputSheet->Cell($outputRow,0) ) < $tol )
		{
		    # Found the right value of x so insert y here
			$outputSheet->Cell($outputRow,$dataColumn) = $yData->[$inputRow];
			$inputRow++;
			$outputRow++;
		}
		elsif ( $inputx < $outputSheet->Cell($outputRow,0) )
		{
		    # This input point is outside the range in the sheet, so skip it and go on to the next
			$inputRow++;
		}
		else
		{
			# Haven't yet come to the right place in the output sheet for this point
			$outputRow++;
		}
	}
}
################################################################################
sub GetPotentialFromForcefield
{
	my ($ff, $interactionType, $sequence, $xVariableName ) =  @_;
	my @xdata = ();
	my @ydata = ();
	my $interactions = $ff->ForcefieldInteractions($interactionType);
	foreach my $ffTerm (@$interactions)
	{
	   if ($ffTerm->Sequence eq $sequence)
	   {
		   foreach my $row (0..$ffTerm->Count-1)
		   {
			   push @xdata, $ffTerm->[$row]->{$xVariableName} ;
			   push @ydata, $ffTerm->[$row]->{"E"} ;
		   }
	   }
	}
	return [\@xdata, \@ydata];
}
################################################################################
sub GetSheet
{
	my ($std, $sheetName) = @_;
	my $sheet;
	eval
	{
		$sheet = $std->Sheets($sheetName);
	};
	return $sheet;
}
################################################################################
sub AddAnalysisResultToSummarySheet
{
	my ($std, $type, $name, $referenceDoc, $chartDoc) = @_;
	
	# revert to master sheet
	$std->ActiveSheetIndex = 0;

	my $row = $std->InsertRow;
	my $tmp = Documents->New($referenceDoc->Name.".xsd");
	$tmp->CopyFrom($referenceDoc);
	KeepSet($tmp,$type);	
	$std->Cell($row,"Structure") = $tmp;
	$tmp->Discard;
	$std->Cell($row,"Name") = $name;
	$std->Cell($row,"Distribution") = $chartDoc;
}	
################################################################################
sub KeepSet
{
	my ($doc,$type) = @_;
	my $sets = $doc->UnitCell->Sets;
	foreach my $set (@$sets) 
	{
		if ($set->Name ne $type)
		{
			$set->Items->Measurements->Delete;
		}
	}
}

################################################################################
sub ColumnHeadingsForIteration
{
	my ($iteration) = @_;
	my $potentialColumnHeading;
	my $pColumnHeading;
	my $pmfColumnHeading;

	if  ( $iteration == 0 )
	{
		$potentialColumnHeading = "Potential";
		$pColumnHeading = "P (target)";
		$pmfColumnHeading = "PMF (target)";
	}
	else
	{
		$potentialColumnHeading = "Potential (trial$iteration)";
		$pColumnHeading = "P (trial$iteration)";
		$pmfColumnHeading = "PMF (trial$iteration)";
	}
	
	return ( $potentialColumnHeading, $pColumnHeading, $pmfColumnHeading );
}
################################################################################
sub EstimatePotential
{
	my ($std,$type,$iteration, $relaxationFactor, $probabilityCutoff) = @_;
	
	# collect the data for fitting
	my $sheet = $std->Sheets($type);
	
	my ( $potentialColumnHeading, $probabilityColumnHeading, $pmfColumnHeading ) = ColumnHeadingsForIteration($iteration);

	my $maxProbability = GetMaximumProbabilityValue($sheet, $probabilityColumnHeading); # based on the current trial trajectory
	my $probabilityThreshold = $probabilityCutoff*$maxProbability;
    
	my @xData = ();
	my @yData = ();
	for(my $row = 0; $row < $sheet->RowCount; $row++)
	{
		my $probabilityTarget = $sheet->Cell($row,"P (target)");
		my $pmfTarget = $sheet->Cell($row,"PMF (target)");
		my $pmfTrial;
		eval { $pmfTrial = $sheet->Cell($row, $pmfColumnHeading )};
		my $potentialTrial;
		eval { $potentialTrial = $sheet->Cell($row, $potentialColumnHeading) };

		my $E;
		if (!defined $probabilityTarget || $probabilityTarget eq "" || $probabilityTarget < $probabilityThreshold )
		{
			# Outside the range that we can estimate, so leave E undefined
		}
		elsif (defined $pmfTarget && ($pmfTarget ne "") && defined $pmfTrial && ($pmfTrial ne "") && defined $potentialTrial)
		{
			# Correct the potential that we used last time using the newly calculated pmf
			$E = $potentialTrial + $relaxationFactor * ( $pmfTarget - $pmfTrial);
			$E = $E * 1.0; # multiply by 1 to force to floating point (avoids a problem with type conversion in the Perl - Materials Studio binding);
		}
		elsif (defined $pmfTarget && ($pmfTarget ne ""))
		{
			# On the first iteration, use the target pmf as the trial potential
			$E = $pmfTarget * 1.0; # multiply by 1 to force to floating point (avoids a problem with type conversion in the Perl - Materials Studio binding);
		}

			push @xData, $sheet->Cell($row,0);		
		push @yData, $E; 
	}
	
	return [\@xData,  \@yData];
}

################################################################################
sub ZeroPotentialAfterNthZero
{
	my ($potential, $numZeroes) = @_;
	
	my $xData = $potential->[0];
	my $yDataIn = $potential->[1];
	my @yDataOut = ();
	my $numZeroesSoFar = 0;
	for my $row (0..$#$xData )
	{
		my $yValue = $yDataIn->[$row];
		# Even number of zeroes so far => potential should be positive
		# Odd number of zeroes so far => potential should be negative
		if ( ($numZeroesSoFar%2 == 0) && $yValue < 0 )
		{
			$numZeroesSoFar = $numZeroesSoFar + 1;
		}
		elsif ( ($numZeroesSoFar%2 == 1) && $yValue > 0 )
		{
			$numZeroesSoFar = $numZeroesSoFar + 1;
		}
		
		if ( $numZeroesSoFar < $numZeroes )
		{
			push @yDataOut, $yValue * 1.0; # multiply by 1 to force to floating point (avoids a problem with type conversion in the Perl - Materials Studio binding)
		}
		else
		{
			push @yDataOut, 0;
		}		
	}
	
	return [ [ @$xData ],  \@yDataOut ];
}

################################################################################
sub GetMaximumProbabilityValue
{
	my ($sheet, $yVariableName) = @_;
	
	my $rowCount = $sheet->RowCount;
	
	# find the maximum probability
	my $maxProb = 0;
	for(my $row = 0; $row < $rowCount; $row++)
	{
		my $Prob = $sheet->Cell($row,$yVariableName);
		$maxProb = $Prob if (defined $Prob && $Prob > $maxProb);
	}
	if ( $maxProb == 0 ) { die "Probability distribution is null\n" };
    return $maxProb;
}
################################################################################
sub FillGapsInPotential
{
	my ($data, $mode ) = @_;
	
	my $isPeriodic = ($mode eq "Periodic");
	
	my $X = $data->[0];
	my $E = $data->[1];
	
	my $N = scalar(@$E);
	if  ($isPeriodic)
		{
		$N = $N - 1; # last point is periodic image of first so ignore for now
		}
	
	my ($lastDefined, $nextDefined) = FindNearestDefinedValuesForEachPoint($data, $isPeriodic);
	
	for my $i (0..$N-1)
	{
		if (!defined $E->[$i])
		{
			my $iLast = $lastDefined->[$i];
			my $iNext = $nextDefined->[$i];
			if ( defined $iLast && defined $iNext )
		{
					# Linear interpolation between the nearest known points
					my $Elast = $E->[ $iLast % $N ];
					my $Enext = $E->[ $iNext % $N ];
					$E->[$i] = ( $Elast * ($iNext-$i) + $Enext * ($i - $iLast) ) / ( $iNext - $iLast );
		}
			elsif ( defined $iLast )
			{
				# Beyond the right hand edge of the known data, so we need to extrapolate
				if ($mode eq "LinearAtLargeX")
				{
					# Use the nearest known point and extrapolate linearly
					$E->[$i] = $E->[ $iLast % $N ] * ( 1 + ($i-$iLast)/$N );
	}
				elsif ($mode eq "QuadraticAtEnds")
				{
					# Extrapolate quadratically from the estimated middle position of the known data
					# through the last known point
					my $iMidPoint = 0.5 * ($nextDefined->[0] + $iLast);
					my $width = $iLast-$iMidPoint;
					if ($width <= 0)
					{
						$width = 1;
					}
					$E->[$i] = $E->[ $iLast % $N ] * (($i-$iMidPoint) / $width)**2;
				}
				elsif ($mode eq "TruncateAtLargeX")
	{
					# Don't need to do anything
				}
				else
		{
							# Default behaviour
							# Copy the value from the nearest known point
							$E->[$i] = $E->[ $iLast % $N ];
				}
		}
			elsif ( defined $iNext )
		{
				# Beyond the left hand edge of the known data, so we need to extrapolate
				if ($mode eq "QuadraticAtEnds")
				{
					# Extrapolate quadratically from the estimated middle position of the known data
					# through the first known point
					my $iMidPoint = 0.5 * ($lastDefined->[-1] + $iNext);
					my $width = $iMidPoint - $iNext;
					if ($width <= 0)
					{
						$width = 1;
					}
					$E->[$i] = $E->[ $iNext % $N ] * ( ($iMidPoint-$i) / $width)**2;
				}
				else
				{
				# Copy the value from the nearest known point
				$E->[$i] = $E->[ $iNext % $N ];
			}
		}
	}
	}
	
	if ($isPeriodic)
	{
	    # Force periodicity if required
		$E->[$N] = $E->[0];
	}
}
################################################################################
sub FindNearestDefinedValuesForEachPoint
{
	my ($data, $isPeriodic) = @_;
	
	my $E = $data->[1];
	my $N = scalar(@$E);
	
	# Pass through the sequence forwards to find the nearest defined value in the backward direction
	my $startF = 0;
	if  ($isPeriodic)
	{
		$N = $N - 1; # last point is periodic image of first so ignore
		$startF = -$N+1;
	}
	my $endF = $N;

	my @lastDefined;
	my $lastDefined;
	for my $i ($startF..$endF)
	{
		if (defined $E->[$i % $N])
		{
			$lastDefined = $i;
		}
		if ($i >= 0 && $i < $N)
	{
			$lastDefined[$i] = $lastDefined;
		}
	}

	# Pass through the sequence backwards to find the nearest defined value in the forward direction
	my $endB = $N-1;
	if  ($isPeriodic)
	{
		$endB = 2*$N-1;
	}
	my $startB = 0;
	my $nextDefined;
	my @nextDefined;

	for my $i (reverse($startB..$endB))
	{
		if (defined $E->[$i % $N])
		{
			$nextDefined = $i;
		}
		if ($i >= 0 && $i < $N)
	{
			$nextDefined[$i] = $nextDefined;
		}
	}
	
	return (\@lastDefined, \@nextDefined);
}
