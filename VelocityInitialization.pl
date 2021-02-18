#!perl

use strict;
use MaterialsScript qw(:all);

# This script contains a number of functions which together provide the implementation
# of a function InitializeAtomVelocities, which set velocities on atoms according to
# a specified temperature.
#
# This can be used to provide alternative starting configurations for a molecular dynamics
# calculation. For example, velocities could be set according to temperature, then those
# of a particular molecule could be given an impulse. Alternatively, two different
# temperatures could be used for different regions, to investigate thermal conductivity.

sub RandomGaussian {
	# The routine returns numbers randomised within the normal distribution with mean 0.0, SD 1.0.
	# It uses the Box-Muller algorithm to generate these from a uniform distribution.
	my $x = 0.0;
	my $y = 0.0;
	my $sumsq = 0.0;
	while ( $sumsq <= 0.0 || $sumsq > 1.0 )  {
		$x = rand() * 2.0 - 1.0;
		$y = rand() * 2.0 - 1.0;
		$sumsq = $x * $x + $y * $y ;
	}
	my $factor = sqrt ( -2.0 * log($sumsq) / $sumsq );
	my $retval = $x * $factor;
	return $retval;
}

sub BoltzmannFactor {
	# returns the Boltzmann constant, adjusted to take account of the units that we are
	# working in:
	# mass in amu
	# velocity in Angstrom/picosecond
	my $Boltzmann = 1.380658e-23;
	my $amu = 1.6605402e-27;
	my $velFactor = 1.0e2;
	my $factor = $Boltzmann / ( $amu * $velFactor * $velFactor );
	
	return $factor;
}

sub RandomizeVelocitiesToTemperature {
	# This randomizes the velocities on all the supplied atoms according to a
	# temperature-based Gaussian distribution.

	my ( $atoms, $temperature ) = @_;

	my $BoltzMann = BoltzmannFactor ( );

	foreach my $atom (@$atoms) {
		my $mass = $atom->Mass;
		my $tempFactor = sqrt ( $BoltzMann * $temperature / $mass );
		my $velx = RandomGaussian() * $tempFactor;
		my $vely = RandomGaussian() * $tempFactor;
		my $velz = RandomGaussian() * $tempFactor;
		$atom->Velocity = Point(X => $velx, Y => $vely, Z => $velz);
	}
}

sub RemoveLinearMomentum {
	# Removes any net linear momentum from the supplied atoms

	my ( $atoms) = @_;

	# calculate net linear momentum
	
	my $sumx = 0.0;
	my $sumy = 0.0;
	my $sumz = 0.0;
	my $totalMass = 0.0;
	foreach my $atom (@$atoms) {
		my $mass = $atom->Mass;
		$sumx = $sumx + $mass * $atom->Velocity->X;
		$sumy = $sumy + $mass * $atom->Velocity->Y;
		$sumz = $sumz + $mass * $atom->Velocity->Z;
		$totalMass = $totalMass + $mass;
	}

	# convert to average
	$sumx = $sumx / $totalMass; 
	$sumy = $sumy / $totalMass;
	$sumz = $sumz / $totalMass; 

	# and subtract from the atom velocities
	foreach my $atom (@$atoms) {
		my $velx = $atom->Velocity->X - $sumx; 
		my $vely = $atom->Velocity->Y - $sumy; 
		my $velz = $atom->Velocity->Z - $sumz; 
		$atom->Velocity = Point(X => $velx, Y => $vely, Z => $velz);
	}

}

sub CurrentTemperature {
	# Calculates the current temperature of a crystal defined by a set of atoms

	# This function assumes that the velocity property has been defined on each atom.
	# In the context of this script this is always the case. In the more general case
	# it would be safer to test whether the velocity is defined, along the lines of
	#	my $vel = $atom->Velocity;
	#	if $vel { ...
	#

	my ( $atoms) = @_;

	my $ke = 0.0;
	foreach my $atom (@$atoms) {
		my $mass = $atom->Mass;
		my $vx = $atom->Velocity->X;
		my $vy = $atom->Velocity->Y;
		my $vz = $atom->Velocity->Z;
		$ke = $ke + 0.5 * $mass * ( $vx * $vx + $vy * $vy + $vz * $vz );
	}
	
	# The calculation of degrees of freedom assumes that the structure is periodic
	# and that all atoms are moveable.
	my $dof = 3 * $atoms->Count - 3;

	my $temperature = 2.0 * $ke / ( $dof * BoltzmannFactor() );
	
	return $temperature;
}


sub ScaleVelocitiesToTemperature {
	# Scales the velocities on the supplied atoms so that their temperature is as requested

	my ( $atoms, $temperature ) = @_;

	my $currentTemp = CurrentTemperature ( $atoms );
	my $velocityRatio = sqrt ( $temperature / $currentTemp );
	foreach my $atom (@$atoms) {
		my $velx = $atom->Velocity->X * $velocityRatio; 
		my $vely = $atom->Velocity->Y * $velocityRatio; 
		my $velz = $atom->Velocity->Z * $velocityRatio; 
		$atom->Velocity = Point(X => $velx, Y => $vely, Z => $velz);
	}
	
}

sub InitializeAtomVelocities {
	# This function sets the velocity property on atoms of a P-1 crystal based on 
	# a specified input temperature.

	my ( $doc, $temperature ) = @_;

	# The following call assumes that we are working with a system in which all atoms
	# are moveable. If there are any fixed atoms then the $atoms variable should be
	# set accordingly.

	my $atoms = $doc->UnitCell->Atoms;

	# Initialize the random number seed (randomly)
	srand;

	# Randomly set the velocities for all the atoms from a temperature-derived distribution.

	RandomizeVelocitiesToTemperature ( $atoms, $temperature );

	# Because this is a random process, the resultant structure will not necessarily
	# have an overall temperature as requested. Also there may be some net linear momentum
	# which needs to be removed first.

	RemoveLinearMomentum ( $atoms );

	ScaleVelocitiesToTemperature ( $atoms, $temperature );
	
}



my $ureaDoc = $Documents{"urea-P1.xsd"};


InitializeAtomVelocities ( $ureaDoc, 300 );
