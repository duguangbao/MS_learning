#!perl

use strict;
use MaterialsScript qw(:all);

# Modules: Materials Visualizer, DMol3

# Performs a potential energy surface (PES) scan by varying the named angle and distance 
# monitors between set values. A DMol3 calculation is performed and the results are 
# written out to a study table.

# The input document must have a distance and angle monitor defined with the appropriate
# names entered into the user editable settings.

# The user should also define these monitors as fixed using the Modify Constraints dialog.
# prior to running the script.

##########################################################################################
# User editable settings

my $filename = "water";		# document name for the scan PES
my $quality = "Coarse";		# Calculation quality - should be set to Fine for good results

# Settings for the distance monitor

my $distanceName = "O-H";
my $minDistance = 0.95;
my $maxDistance = 1.10;
my $distanceStep = 0.05;

# Settings for the angle monitor

my $angleName = "H-O-H";
my $minAngle = 80;
my $maxAngle = 90.0;
my $angleStep = 5;

# End user editable settings
##########################################################################################

# Load the document
my $doc = $Documents{"$filename.xsd"};

# Create the study table and column headings
my $studyTable = Documents->New("$filename"."PES.std")->ActiveSheet;
$studyTable->ColumnHeading(0) = "Structure";
$studyTable->ColumnHeading(1) = "Angle";
$studyTable->ColumnHeading(2) = "Distance";
$studyTable->ColumnHeading(3) = "Energy";

# Specify the module and some settings
my $dmol3 = Modules->DMol3;
$dmol3->ChangeSettings([	UseSymmetry	=>	0,
				Quality		=>	$quality
			]);

# Get the monitors using the collection API
my $distanceMonitor = $doc->Distances("$distanceName");
my $angleMonitor = $doc->Angles("$angleName");

# Row counter for the study table so that structures are put in the correct rows.

my $rowCounter = 0 ;

# Main loop which iterates through all the variations of angle and distance, calculating
# the energy each time and populating the study table.

for (my $angle = $minAngle; $angle <= $maxAngle; $angle+=$angleStep) {
	
	# Set the angle to the new value
	$angleMonitor->Angle = $angle;
	
	for (my $distance = $minDistance; $distance <= $maxDistance; $distance+=$distanceStep) {
		
		# Set the distance monitor to the new value
		$distanceMonitor->Distance = $distance;

		print "Processing frame $rowCounter, angle is $angle and distance is $distance\n";

		# Calculate the DMol3 Geometry Optimization task within an eval statement in case it fails
		# Putting it in the eval loop means that if it does fail, it will proceed onto the next
		# frame without stopping the script.
		my $output;
		eval {$output = $dmol3->GeometryOptimization->Run($doc);};
		
		if ($@) {
			# Something bad has happened, move onto the next one
			print "!!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!\n";
			print $@;
		
		} else {
			# Calculation has completed successfully, carry on
			my $energy=$output->TotalEnergy;

			# Put the data in the study table
			$studyTable->Cell($rowCounter, 0) = $doc;
			$studyTable->Cell($rowCounter, 1) = $angle;
			$studyTable->Cell($rowCounter, 2) = $distance;
			$studyTable->Cell($rowCounter, 3) = $energy;
		}
		++$rowCounter;
	}
}

print "The script has completed.\n";

