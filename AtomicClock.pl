#!perl
#
# Purpose: Demonstrate building an arbitrary atomic structure that tells the
#          time. Hit Escape to terminate!

use strict;
use warnings;
use MaterialsScript qw(:all);
use Time::HiRes 'sleep';

#Constants
use constant seconddist   => 4.25;            # Length of second hand
use constant minutedist   => 3.5;             # Length of minute hand
use constant hourdist     => 2.75;            # Length of hour hand
use constant tickdist     => 4.00;            # Distance to center of perimeter ticks
use constant TWOPI        => 2 * 3.141592654; # Number of radians in a circle
use constant BALLANDSTICK => "Ball and Stick";

#Global variables
my $doc        = undef; # The document
my $secondHand = undef; # The second-hand atom in the document
my $minuteHand = undef; # The minute-hand atom in the document
my $hourHand   = undef; # The hour-hand atom in the document

#Calculate the XYZ for a given point on the clock
sub ClockPoint {
    my ($dist, $divs, $val, $z) = @_;
    my $angle = $val / $divs * TWOPI;
    return Point(X => $dist * sin($angle),
                 Y => $dist * cos($angle),
                 Z => $z);
}

#Check if a document exists, return 1 (true) if so 0 (false) otherwise
sub ClockExists {
    eval q($Documents{"AtomicClock.xsd"}); return 0 if $@;
    
    return 1;
}
	
#Create the clock face
sub CreateClock {
    $doc = Documents->New("AtomicClock.xsd");
    
    # Create an atom at the center of the clock
    my $center = $doc->CreateAtom("C", Point(Z => -0.5));
    
    # Create the atoms for the hour markers at the clock perimeter
    for (my $hours=0; $hours<12; ++$hours) {
        $doc->CreateAtom("Si", ClockPoint(tickdist, 12, $hours, -0.5))
            ->Style = BALLANDSTICK;
    }
    
    #Create atoms and bonds for the ends of the second, minute and hour hands
    $secondHand = $doc->CreateAtom("H", Point(Z => +0.5), ["IsVisible" => "No", "Name" => "Second Hand"]);
    $minuteHand = $doc->CreateAtom("H", Point(Z => +0.25),["Name" => "Minute Hand"]);
    $hourHand   = $doc->CreateAtom("O", Point(Z => 0),    ["Name" => "Hour Hand"]);
    $doc->CreateBestFitLine([$center, $secondHand],       ["LineExtentPadding" => 0]);
    $doc->CreateBond($center, $minuteHand, "Single")->Style = BALLANDSTICK;
    $doc->CreateBond($center, $hourHand,   "Double")->Style = BALLANDSTICK;
    
    $minuteHand->Style = BALLANDSTICK;
    $hourHand->  Style = BALLANDSTICK;
}

#Map variables required to update the existing clock
sub AcquireClock {
    $doc = $Documents{"AtomicClock.xsd"};
    
    $secondHand = $doc->Atoms("Second Hand");
    $minuteHand = $doc->Atoms("Minute Hand");
    $hourHand   = $doc->Atoms("Hour Hand");
}

#Update the position of the clock hands 
sub UpdateClockPosition {
    my @timedata = localtime(time);
    my $seconds  = $timedata[0] + ($timedata[1] * 60) + ($timedata[2] * 3600);
    my $hours    = $seconds / 3600;
    my $mins     = $seconds / 60;
    
    $secondHand->XYZ = ClockPoint(seconddist, 60, $seconds, $secondHand->Z);
    $minuteHand->XYZ = ClockPoint(minutedist, 60, $mins,    $minuteHand->Z);
    $hourHand  ->XYZ = ClockPoint(hourdist,   12, $hours,   $hourHand  ->Z);
}

#main
if (ClockExists) {
    AcquireClock;
} else {
    CreateClock;
}

while (1) {
    UpdateClockPosition;
    $doc->UpdateViews;
    sleep 0.1;
}
