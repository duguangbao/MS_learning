#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);
use Math::MatrixReal;

# Calculates a stress-strain curve using the Souza-Martins barostat

# Applies a stress in one direction (zero stresses in the others) and measures the strain.
# Output is a study table containing the stress-strain curve.
# Intended to be run from the User menu and operates on the Active Document, which should be a 
# 3D periodic structure. 

# Input parameters are:
#		NumZeroStressEquilibrationCycles 零应力平衡循环次数
#		NumTimeStepsPerZeroStressEquilibrationCycle 每个零应力平衡循环的时间步数
#		StressComponent 应力分量
#		Stresses 
#		NumEquilibrationStepsPerStressValue 每个应力值的平衡步数
#		NumProductionStepsPerStressValue 每个应力值的生产步骤数

# Note that there is an initial equilibration phase at zero stress. This is done as a sequence of several cycles, with the cell size and shape
# reset at the beginning of each cycle (to the average size over the previous cycle). (This is intended to accelerate the equilibration by stopping
# oscillations of the cell.) This is controlled by parameters NumZeroStressEquilibrationCycles and NumTimeStepsPerZeroStressEquilibrationCycle.

# The StressComponent is one of XX YY ZZ XY XZ YZ. All other stress components are set to zero. If XX YY or ZZ, then the results corresponds to a
# Young's modulus; if XY, XZ, or YZ then a shear modulus

# Stresses is a list of numbers. This is the sequence of stresses to be applied in units of GPa. At each stress, there is a period of equilibration
# followed by a period of production. This is controlled by parameters NumEquilibrationStepsPerStressValue and NumProductionStepsPerStressValue


# Input parameters
my $inputStructure;
my $numZeroStressEquilibrationCycles; 
my $numTimeStepsPerZeroStressEquilibrationCycle;
my $stressComponent;
my @stresses;
my $numEquilibrationStepsPerStressValue;
my $numProductionStepsPerStressValue;

# This block acquires the input structure and parameters from the User dialog. If running the script directly instead of
# via the User menu, input parameters can be set directly in the following block
eval {

	$inputStructure = Documents->ActiveDocument; 
	
	my %Args;
	GetOptions(\%Args, 	"NumZeroStressEquilibrationCycles=i", 
				"NumTimeStepsPerZeroStressEquilibrationCycle=i",
				"StressComponent=s",
				"Stresses=s",
				"NumEquilibrationStepsPerStressValue=i", 
				"NumProductionStepsPerStressValue=i"
		);
	
	# Initial equilibration phase at zero stress
	$numZeroStressEquilibrationCycles = $Args{"NumZeroStressEquilibrationCycles"};
	$numTimeStepsPerZeroStressEquilibrationCycle = $Args{"NumTimeStepsPerZeroStressEquilibrationCycle"}; 
	
	# Direction of stress
	$stressComponent =  $Args{"StressComponent"}; # One of XX YY ZZ XY XZ YZ
	
	# These stress values are applied in turn and appear in the output
	my $stressString = $Args{"Stresses"};
	@stresses = split( / / , $stressString ); # stresses in GPa.
	
	# Number of time steps at each stress value
	$numEquilibrationStepsPerStressValue = $Args{"NumEquilibrationStepsPerStressValue"}; # Equilibration time steps for each stress 
	$numProductionStepsPerStressValue =    $Args{"NumProductionStepsPerStressValue"};    # Production time steps for each stress
};



# This block sets the input parameters if not running from the User menu
if (!defined $inputStructure)
{
    # These values were chosen for testing purposes. Modify as required
	$inputStructure = $Documents{"Fe.xsd"};
	$numZeroStressEquilibrationCycles = 5;
	$numTimeStepsPerZeroStressEquilibrationCycle = 100;
	$stressComponent = "XX";
	@stresses = (1, 2, 3);
	$numEquilibrationStepsPerStressValue = 500;
	$numProductionStepsPerStressValue = 500;
}

# Stem of the name to be used for output files. Modify if required
my $outputName = $inputStructure->Name;
	
# Energy settings. Modify as required
my $energySettings = [
	CurrentForcefield => "COMPASS"
];

###################

# Columns in the output study table
use constant INPUT_STRESS => 0;
use constant MEAN_STRESS => 1;
use constant ERROR_STRESS => 2;
use constant STRAIN => 3;
use constant MODULUS => 4;

my $std = Documents->New($outputName ."Stress.std");
$std->ColumnHeading(INPUT_STRESS) = "Input stress / GPa";
$std->ColumnHeading(MEAN_STRESS) = "Measured stress / GPa";
$std->ColumnHeading(ERROR_STRESS) = "Error in stress / GPa";
$std->ColumnHeading(STRAIN) = "Strain";
$std->ColumnHeading(MODULUS) = "Modulus / GPa";

my $stressMode = "Cauchy" ;  # Constant stress mode, either "Cauchy" or "2nd Piola-Kirchhoff"
my @cellVectorComponents = qw{ AX AY AZ BX BY BZ CX CY CZ};
my @propertyNames = qw{ AX AY AZ BX BY BZ CX CY CZ StressXX StressYY StressZZ StressXY StressXZ StressYZ };
my $targetTrajectoryFrequency = 1000;

# Equilibrate at zero stress. Do $numZeroStressEquilibrationCycles cycles each of $numTimeStepsPerZeroStressEquilibrationCycle steps
my $startStructure = $inputStructure;
SetCellMatrixVelocityToZero($startStructure);
for my $cycle (1..$numZeroStressEquilibrationCycles)
{
	my $stress = 0;

	my $trajectoryFrequency = EstimateTrajectoryFrequency( $targetTrajectoryFrequency, $numTimeStepsPerZeroStressEquilibrationCycle );
	my $results = RunSimulation( $startStructure, [ Stress => $stress,
							       ReferenceCell => undef,
							       OutputName => $outputName . "Equil",
							       TrajectoryFrequency => $trajectoryFrequency,
							       NumberOfSteps => $numTimeStepsPerZeroStressEquilibrationCycle ]);
								   
	my ($averages) = CalculateAverages($results, 0);
	my $outputStructure = $results->Structure;
	
	# Copy cell parameters from $averages to %averageCell
	my $averageCell = GetCellMatrixFromAveragedData( $averages );
	
	$startStructure = $outputStructure;
	
	# The starting cell matrix on the next cycle is the average from this one
	SetCellMatrix( $startStructure, $averageCell );
	SetCellMatrixVelocityToZero($startStructure);
}




# Now do a run at zero stress (to get the cell dimensions)
my $stress = 0;
my $initialCell = GetCellMatrix($startStructure);
my $trajectoryFrequency = EstimateTrajectoryFrequency( $targetTrajectoryFrequency, $numProductionStepsPerStressValue );
my $results = RunSimulation( $startStructure, [ Stress => $stress,
							   ReferenceCell => $initialCell,
							   OutputName => $outputName,
							   TrajectoryFrequency => $trajectoryFrequency,
							   NumberOfSteps => ($numEquilibrationStepsPerStressValue + $numProductionStepsPerStressValue) ]);
								  
my ($averages) = CalculateAverages($results, $numEquilibrationStepsPerStressValue/$trajectoryFrequency);

# Copy cell parameters from $averages to %referenceCell
my $referenceCell = GetCellMatrixFromAveragedData($averages);

my $strainTensor = CalculateSymmetricStrainTensor( $referenceCell, $referenceCell );
my $outputStructureZeroStress = $results->Structure;
AddToStudyTable($std, $averages, 0, 0, $strainTensor);


# Now do non-zero stresses. Note that for each stress value, we do a single dynamics run to do both equilibration and production.
# The equilibration period will be ignored by the analysis
$startStructure = $outputStructureZeroStress;
for my $point (0..scalar(@stresses)-1)
{
	my $stress = $stresses[$point];
	my $trajectoryFrequency = EstimateTrajectoryFrequency( $targetTrajectoryFrequency, $numProductionStepsPerStressValue );
	my $results = RunSimulation($startStructure, [ Stress => $stress,
								  ReferenceCell => $referenceCell,
								  OutputName => $outputName,
								  TrajectoryFrequency => $trajectoryFrequency,
								  NumberOfSteps => ($numEquilibrationStepsPerStressValue + $numProductionStepsPerStressValue) ]);
	
	my ($averages) = CalculateAverages($results, $numEquilibrationStepsPerStressValue/$trajectoryFrequency);
	my $averageCell = GetCellMatrixFromAveragedData($averages);
	my $strainTensor = CalculateSymmetricStrainTensor( $referenceCell, $averageCell );
	my $outputStructure = $results->Structure;	
	AddToStudyTable($std, $averages, $point+1, $stress, $strainTensor);
	$startStructure = $outputStructure;
}


####

# Adds a row to the study table for one value of applied stress.
# Input is: the study table; a hash containing averages (and standard errors) of the measured stress components;
#           row number; the applied stress; the average strain tensor (as a hash)
sub AddToStudyTable
{
	my ($std, $averages, $row, $stress, $strainTensor) = @_;
	my $meanStress = $averages->{"MeanStress" . $stressComponent}; 
	$std->Cell($row,INPUT_STRESS) = $stress + 0; # +0 to force string to number
	$std->Cell($row,MEAN_STRESS) = $meanStress; 
	$std->Cell($row,ERROR_STRESS) = $averages->{"ErrorStress" . $stressComponent};       

	# Calculate strain and modulus
	my $tensorialStrain = $strainTensor->{$stressComponent};
	my $engineeringStrain;
	if ( $stressComponent eq "XY" or $stressComponent eq "XZ" or $stressComponent eq "YZ" )
	{
		$engineeringStrain = 2 * $tensorialStrain;
	}
	elsif ( $stressComponent eq "XX" or $stressComponent eq "YY" or $stressComponent eq "ZZ" )
	{
		$engineeringStrain = $tensorialStrain;
	}
	
	$std->Cell($row,STRAIN) = $engineeringStrain;
	if ($engineeringStrain != 0)
	{
		$std->Cell($row,MODULUS) = $meanStress / $engineeringStrain;
	}
}


# Does a dynamics run. Input is a structure document; a hash containing parameter values for
# Stress, ReferenceCell, OutputName, TrajectoryFrequency, NumberOfSteps
# ReferenceCell is supplied as a hash { AX => ??, AY => ??, ...,  CZ => ?? }
sub RunSimulation
{
	my ($startStructure,$parameterList) = @_;
	my %parameters = @$parameterList;
	my $stress = $parameters{"Stress"};
	my $referenceCell = $parameters{"ReferenceCell"};
	my $outputName = $parameters{"OutputName"};
	my $trajFreq = $parameters{"TrajectoryFrequency"};
	my $numSteps = $parameters{"NumberOfSteps"};
	
	my $filename = $outputName . $stress;
	$filename =~ s/\./p/; # Avoid filenames with more than one '.' character
	
	my $doc = Documents->New( $filename . ".xsd");
	$doc->CopyFrom( $startStructure );

	Modules->Forcite->ChangeSettings([
		Ensemble3D => "NPT",
		Barostat => "Souza-Martins",
		SouzaMartinsBarostatCellTimeConstant => 10,
		SouzaMartinsConstantStressMode => $stressMode,
		Thermostat => "NHL",
		TrajectoryFrequency => $trajFreq,
		NumberOfSteps => $numSteps,
		"Stress" . $stressComponent => $stress,
		TimeStep => 1,
		WriteLevel => "Silent"
	]);

	Modules->Forcite->ChangeSettings( $energySettings );

	# Extracts the reference cell components AX...CZ from the hash and sets them into Forcite as properties SouzaMartinsReferenceCellVectorAX etc
	# Note that these are relevant only for 2nd Piola-Kirchhoff stress mode
	if (defined $referenceCell)
	{
		map { Modules->Forcite->ChangeSettings([ "SouzaMartinsReferenceCellVector" . $_ => $referenceCell->{$_} ]) } @cellVectorComponents;
	}
	
	my $results = Modules->Forcite->Dynamics->Run($doc);
	return $results;
}


# Takes a dynamics result object and computes means and standard errors of the cell matrix and stress components.
# The first  $numEquilibrationFrames of the trajectory are discarded.
# Output is a hash { MeanAX => ???, ErrorAX => ???, ... }
sub CalculateAverages
{
	my ($results, $numEquilibrationFrames) = @_;
	my $traj = $results->Trajectory;
	my $numFrames = $traj->NumFrames - $numEquilibrationFrames;

	if ($numFrames <= 0 )
	{
		die ("No production frames in trajectory\n");
	}
	
	# Ignore the first $numEquilibrationFrames of the trajectory, then divide the rest into blocks.
	# We use the standard deviation of the block averages to estimate standard error
	
	# First compute the total and total square of the block averages (of the stresses and cell matrix components)
	my $numBlocks = 4;
	my %averages;
	for my $block (0..($numBlocks-1))
	{
		my $startFrame = int($block * $numFrames / $numBlocks) + $numEquilibrationFrames;
		my $nextStartFrame = int(($block + 1) * $numFrames / $numBlocks) + $numEquilibrationFrames;
		my $data = AverageStress($startFrame, $nextStartFrame-$startFrame, $traj  ); # $data is a hash containing block averages
		for my $property (@propertyNames)
		{
			$averages{"Total" . $property} += $data->{$property};
			$averages{"TotalSqd" . $property} += $data->{$property}**2;
		}	
	}

	# Finally, the mean and standard error 
	for my $property (@propertyNames)
	{
		my $mean = $averages{"Total" . $property} / $numBlocks;
		my $variance = $averages{"TotalSqd" . $property}/$numBlocks - $mean**2;
		
		$averages{"Mean" . $property} = $mean;
		$averages{"Error" . $property} = sqrt( $variance * $numBlocks / ($numBlocks - 1));
	}

	return (\%averages);
}

# Takes a trajectory and returns the average of the stress tensor (and the cell matrix)
# over some of its frames. Result is a hash { StressXX => ???, ... , AX => ??? ... }
sub AverageStress
{
	my ($startFrame, $numFrames, $traj) = @_;

	my %cumulativeData;	
	for my $frame ($startFrame..($startFrame+$numFrames-1))
	{
		$traj->CurrentFrame = $frame;

		my $stressTensor = GetStressTensor($traj);
		for my $key (keys %$stressTensor)
		{
			$cumulativeData{$key} += $stressTensor->{$key};
		}
		
		my $cellMatrix = GetCellMatrix($traj);
		for my $key (keys %$cellMatrix)
		{
			$cumulativeData{$key} += $cellMatrix->{$key};
		}		
	}

	my %results;
	for my $key (keys %cumulativeData)
	{
		$results{$key} = $cumulativeData{$key} / $numFrames;
	}	

	return \%results;
}

# Gets the stress tensor from a structure and returns it as a hash { StressXX => ?? , StressXY => ?? , ... } 
sub GetStressTensor
{
	my ($structure) = @_;
	my $symmetrysystem = $structure->SymmetrySystem;
	my $stress = $symmetrysystem->Stress;

	my %stressTensor;
	$stressTensor{"StressXX"} = $stress->Eij(1,1);
	$stressTensor{"StressYY"} = $stress->Eij(2,2);
	$stressTensor{"StressZZ"} = $stress->Eij(3,3);
	$stressTensor{"StressXY"} = $stress->Eij(1,2);
	$stressTensor{"StressXZ"} = $stress->Eij(1,3);
	$stressTensor{"StressYZ"} = $stress->Eij(2,3);		

	return \%stressTensor;
}


# Gets the cell matrix from a structure and returns it as a hash { AX => ?? , AY => ?? , ... CZ => ?? } 
sub GetCellMatrix
{
	my ($structure) = @_;
	my $symmetrysystem = $structure->SymmetrySystem;
	my $symmdef = $symmetrysystem->SymmetryDefinition;
	my %cellMatrix;
	$cellMatrix{"AX"} = $symmdef->VectorA->X;
	$cellMatrix{"AY"} = $symmdef->VectorA->Y;
	$cellMatrix{"AZ"} = $symmdef->VectorA->Z;
	$cellMatrix{"BX"} = $symmdef->VectorB->X;
	$cellMatrix{"BY"} = $symmdef->VectorB->Y;
	$cellMatrix{"BZ"} = $symmdef->VectorB->Z;
	$cellMatrix{"CX"} = $symmdef->VectorC->X;
	$cellMatrix{"CY"} = $symmdef->VectorC->Y;
	$cellMatrix{"CZ"} = $symmdef->VectorC->Z;
	
	return \%cellMatrix;
}



# Sets the cell matrix of a structure from a hash { AX => ?? , AY => ?? , ... CZ => ?? } 
sub SetCellMatrix
{
	my ($structure, $cellMatrix) = @_;
	my $symmetrysystem = $structure->SymmetrySystem;
	my $symmdef = $symmetrysystem->SymmetryDefinition;

	$symmdef->VectorA->X = $cellMatrix->{"AX"};
	$symmdef->VectorA->Y = $cellMatrix->{"AY"};
	$symmdef->VectorA->Z = $cellMatrix->{"AZ"}; 
	$symmdef->VectorB->X = $cellMatrix->{"BX"}; 
	$symmdef->VectorB->Y = $cellMatrix->{"BY"}; 
	$symmdef->VectorB->Z = $cellMatrix->{"BZ"}; 
	$symmdef->VectorC->X = $cellMatrix->{"CX"}; 
	$symmdef->VectorC->Y = $cellMatrix->{"CY"}; 
	$symmdef->VectorC->Z = $cellMatrix->{"CZ"}; 
}

# Sets the cell matrix rate of change of a structure to zero
sub SetCellMatrixVelocityToZero
{	
	my ($structure) = @_;
	my $symmetrysystem = $structure->SymmetrySystem;
	my $symmdef = $symmetrysystem->SymmetryDefinition;
	
	$symmdef->VectorAVelocity->X = 0;
	$symmdef->VectorAVelocity->Y = 0;
	$symmdef->VectorAVelocity->Z = 0;
	$symmdef->VectorBVelocity->X = 0;
	$symmdef->VectorBVelocity->Y = 0;
	$symmdef->VectorBVelocity->Z = 0;
	$symmdef->VectorCVelocity->X = 0;
	$symmdef->VectorCVelocity->Y = 0;
	$symmdef->VectorCVelocity->Z = 0;
	
}

# Gets a cell matrix as a hash: { AX => ???, AY => ???, ... , CZ => ??? }
# from the averaged data (a hash of the form { MeanAX => ??? etc } )
sub GetCellMatrixFromAveragedData
{
	my ($averages) = @_;
	my %referenceCell;
	map{ $referenceCell{$_} = $averages->{"Mean" . $_ } }  @cellVectorComponents;
	return \%referenceCell;
}



# Calculates the strain tensor (1/2) * ( h . h0^(-1)  +  h0^(-T) . h^T ) - 1 
# where h is the (strained) cell matrix and h0 is the reference cell matrix.
# Input cell matrices are hashes { AX => ???, AY => ???, ... , CZ => ??? }
# Output matrix is hash          { XX => ???, XY => ???, ... , ZZ => ??? }
sub CalculateSymmetricStrainTensor
{
	my ($referenceCell, $currentCell) = @_;
	my $currentCellAsMatrixReal = CellMatrixHashToMatrixReal($currentCell);
	my $referenceCellAsMatrixReal = CellMatrixHashToMatrixReal($referenceCell);
	my $inverseReferenceCell = $referenceCellAsMatrixReal->inverse();
	my $m1 = $currentCellAsMatrixReal->multiply( $inverseReferenceCell );

	my $identity = new Math::MatrixReal(3,3);
	$identity->one();
	
	my $m2 = new Math::MatrixReal(3,3);
	$m2->subtract($m1,$identity);
	
	my $m2transpose = new Math::MatrixReal(3,3);
	$m2transpose->transpose($m2);
	
	my $m3 = $m2 + $m2transpose;
	my $m4 = new Math::MatrixReal(3,3);
	$m4->multiply_scalar($m3, 0.5);
	return MatrixRealToCartesianHash($m4);
}


# Converts a matrix from a hash { AX => ??? etc } to Math::MatrixReal
sub CellMatrixHashToMatrixReal
{
	my ($matrixAsHash) = @_;
	my $m = Math::MatrixReal->new_from_rows( [ [ $matrixAsHash->{"AX"}, $matrixAsHash->{"BX"}, $matrixAsHash->{"CX"} ],
	                                           [ $matrixAsHash->{"AY"}, $matrixAsHash->{"BY"}, $matrixAsHash->{"CY"} ],
	                                           [ $matrixAsHash->{"AZ"}, $matrixAsHash->{"BZ"}, $matrixAsHash->{"CZ"} ] ] );
	return $m;
}

# Converts a matrix from a Math::MatrixReal to a hash { XX => ??? etc } 
sub MatrixRealToCartesianHash
{
	my ($matrixReal) = @_;
	my %matrixAsHash;
	$matrixAsHash{"XX"} = $matrixReal->element(1,1);
	$matrixAsHash{"XY"} = $matrixReal->element(1,2);
	$matrixAsHash{"XZ"} = $matrixReal->element(1,3);
	$matrixAsHash{"YX"} = $matrixReal->element(2,1);
	$matrixAsHash{"YY"} = $matrixReal->element(2,2);
	$matrixAsHash{"YZ"} = $matrixReal->element(2,3);
	$matrixAsHash{"ZX"} = $matrixReal->element(3,1);
	$matrixAsHash{"ZY"} = $matrixReal->element(3,2);
	$matrixAsHash{"ZZ"} = $matrixReal->element(3,3);
	return \%matrixAsHash;
}

# Gives a sensible value for TrajectoryFrequency, given a suggested $targetFrequency.
# Whichever is smaller out of $targetFrequency and $numberOfSteps / 10
sub EstimateTrajectoryFrequency
{
	my ($targetFrequency, $numberOfSteps) = @_;
	my $trajectoryFrequency = $targetFrequency; 
	if ( $numberOfSteps < 10 * $targetFrequency )
	{
		$trajectoryFrequency = int( $numberOfSteps / 10 );
	}
	if ( $trajectoryFrequency < 1 )
	{
		$trajectoryFrequency = 1;
	}
	return $trajectoryFrequency;

}
