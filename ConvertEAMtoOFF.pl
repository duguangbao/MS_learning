###############################################################################
#!perl
#
# Description: Converts EAM tabulated forcefields to OFF
#
# Author: BIOVIA
# Version: 1.0 (August 2017)
# Required: Materials Studio 2018 or later.
#
# This script converts tabulated embedded atom model (EAM) potentials to 
# .off format, for use in the Forcite module of Materials Studio. 
#
# The following formats are supported:
# 1. eam.fs
# 2. eam.alloy
# The file type is parsed from the length of the file.
#
# Required input:
# - A text document containing the EAM data
# 
# Output documents:
# - A forcefield document containing the converted EAM data
###############################################################################

use strict;
use warnings;
use MaterialsScript qw(:all);
use constant { kcal_per_mol_per_eV => 23.06054721 };
use constant { FS => 0, ALLOY => 1};
{
################################################################################
#  BEGIN USER INPUT                                                            #

	# set the name of text file containing the EAM data
	my $fileName = "CuNi_Example.eam.fs";  # without the .txt extension
	
#  END USER  INPUT                                                             #
################################################################################
	Run($fileName);
}

sub Run
{
	my ($fileName) = @_;
	my $txt;
	
	eval
	{
		$txt = $Documents{"$fileName.txt"};
		print "Opening file: $fileName.txt\n";
	};
	if($@)
	{
		die "The file $fileName.txt specified for input does not exist";
	}

	my $fileFormat = GetFileFormatFromLines($txt);

	# Create a new forcefield document
	my $off = Documents->New("$fileName.off");

	print "Creating new forcefield document: $fileName.off\n";

	# switch on the EAM interactions in the forcefield
	$off->UseEAM = "Yes";

	# use the first 3 lines as description
	my $lines = $txt->Lines;
	my $row = 0;
	$off->Description = $lines->Item($row++) . "\n";
	$off->Description .= $lines->Item($row++) . "\n";
	$off->Description .= $lines->Item($row++);

	# parse the line containing the elements
	my @tokens = split ' ', $lines->Item($row++);
	my $numElements = $tokens[0];
	my @elements;
	for(my $i = 1; $i < scalar(@tokens); $i++)
	{
		my $element = $tokens[$i];
		push @elements, $element;
		$off->CreateType($element, [ElementSymbol => $element]);
	}

	# parse the line containing the density
	@tokens = split ' ', $lines->Item($row++);
	my $numDensity = 1.0 * $tokens[0];
	my $densityBin = 1.0 * $tokens[1];
	my $numRadius = 1.0 * $tokens[2];
	my $radiusBin = 1.0 * $tokens[3];
	my $radiusMax = 1.0 * $tokens[4];

	# run over each element
	for(my $iElementIndex = 0; $iElementIndex < scalar(@elements); $iElementIndex++)
	{
		my $element = $elements[$iElementIndex];
		print "Parsing data for element: $element\n", ;
		
		# parse the line containing the element data
		@tokens = split ' ', $lines->Item($row++);

		my $atomicNumber = $tokens[0];
		my $molecularWeight = $tokens[1];
		my $latticeConstant = $tokens[2];
		my $latticeType = $tokens[3];

		# the embedding functions F(rho)
		$row = ParseEmbeddingFunction($off, $element, $lines, $row, $numDensity, $densityBin);
		#ParseEmbeddingFunction($off, $element, $lines, $row, $numDensity, $densityBin);
		
		# The density functions rho(r). Depending on the fileformat these may be asymmetric (FS, SET) or symmetric (ALLOY)
		
		# Symmetric density functions
		if ($fileFormat eq FS) 
		{
			for(my $jElementIndex = 0; $jElementIndex < scalar(@elements); $jElementIndex++)
			{
				my $element2 = $elements[$jElementIndex];
				$row = ParseElectronDensity($off, $element, $element2, $lines, $row, $numRadius, $radiusBin);
			}
		}
		# Asymmetric density functions
		elsif ($fileFormat eq ALLOY) 
		{
			$row = ParseElectronDensities($off, $element, \@elements, $lines, $row, $numRadius, $radiusBin);		 
		}
	}
	
	# the pair potentials phi(r) (symmetric)
	for(my $iElementIndex = 0; $iElementIndex < scalar(@elements); $iElementIndex++)
	{
		my $element1 = $elements[$iElementIndex];
		for(my $jElementIndex = 0; $jElementIndex <= $iElementIndex; $jElementIndex++)
		{
			my $element2 = $elements[$jElementIndex];
			$row = ParseEAMPairPotential($off, $element1, $element2, $lines, $row, $numRadius, $radiusBin);
		}
	}
	
	$off->Save;
	print "Done.";
}

sub ParseEmbeddingFunction
{
	my ($off, $element, $lines, $start, $numDensity, $densityBin) = @_;
	my $row = $start;
	my @E;
	my @rho;
	my $density = 0;		
	do
	{
		my @tokens = split ' ', $lines->Item($row++);
			
		for(my $j = 0; $j < scalar(@tokens); $j++)
		{
			push @E, $tokens[$j]*kcal_per_mol_per_eV;
			push @rho, $density;
			$density += $densityBin;
		}		
	}
	while(scalar(@rho) < $numDensity);
	
	$off->CreateTerm("Embedding Function",$element,"Tabulated",[rho => @rho, E => @E]);
	
	return $row;
}		

sub ParseElectronDensity
{
	my ($off, $element1, $element2, $lines, $start, $numRadius, $radiusBin) = @_;
	my $row = $start;
	
	my @rho;
	my @R;
	my $radius = 0;
	do
	{
		my @tokens = split ' ', $lines->Item($row++);
			
		for(my $j = 0; $j < scalar(@tokens); $j++)
		{
			next if (scalar(@R) > $numRadius);
			push @rho, $tokens[$j]+0.0;
			push @R, $radius;
			$radius += $radiusBin;
		}		
	}
	while(scalar(@R) < $numRadius);
	
	my $typeSequence = $element1 . "," . $element2;
	$off->CreateTerm("Electron Density",$typeSequence,"Tabulated",[R => @R, rho => @rho]);
	
	return $row;
}

sub ParseElectronDensities
{
	my ($off, $element1, $elements, $lines, $start, $numRadius, $radiusBin) = @_;

	my $row = $start;
	my @elements = @$elements;
	my @rho;
	my @R;
	my $radius = 0;
	do
	{
		my @tokens = split ' ', $lines->Item($row++);
				
		for(my $j = 0; $j < scalar(@tokens); $j++)
		{
			next if (scalar(@R) > $numRadius);
			push @rho, $tokens[$j]+0.0;
			push @R, $radius;
			$radius += $radiusBin;
		}		
	}
	while(scalar(@R) < $numRadius);
	
	for(my $jElementIndex = 0; $jElementIndex < scalar(@elements); $jElementIndex++) 
	{	
		my $typeSequence = $elements[$jElementIndex] . "," . $element1;
		$off->CreateTerm("Electron Density",$typeSequence,"Tabulated",[R => @R, rho => @rho]);
	}
	
	return $row;
}


sub ParseEAMPairPotential
{
	my ($off, $element1, $element2, $lines, $start, $numRadius, $radiusBin) = @_;
	my $row = $start;
	my @E;
	my @R;
	my $radius = 0;
	my $skip = 0;
	do
	{
		my @tokens = split ' ', $lines->Item($row++);
			
		for(my $j = 0; $j < scalar(@tokens); $j++)
		{
			# Note: what is tabulated is phi(r)*r, hence divide by r to get phi(r)
			# To avoid division by zero for r = 0, skip this data point
			if($radius > 0)
			{
				push @E, $tokens[$j]*kcal_per_mol_per_eV/$radius;
				push @R, $radius;
			}
			else
			{
				$skip++;
			}
			
			$radius += $radiusBin;
		}		
	}
	while(scalar(@R) < $numRadius-$skip);
	
	my $typeSequence = $element1 . "," . $element2;		
	$off->CreateTerm("EAM Pair Potential",$typeSequence,"Tabulated",[R => @R, E => @E]);	
	
	return $row;		
}

sub GetFileFormatFromLines
{
	my ($txt) = @_;
	
	print("Checking format: \n");
		
	my $lines = $txt->Lines;
	my $row = 5;
	# parse the line containing the elements
	my @tokens = split ' ', $lines->Item(3);
	my $numElements = $tokens[0];
	my @elements;
	for(my $i = 1; $i < scalar(@tokens); $i++)
	{
		my $element = $tokens[$i];
		push @elements, $element;
	}

 	my @atomicNumbers;
	for (my $i = 0; $i < scalar(@elements); $i++) {
		$atomicNumbers[$i] = AtomicNumber($elements[$i]); 
	}

	# parse the line containing the density
	@tokens = split ' ', $lines->Item(4);
	my $numDensity = 1.0 * $tokens[0];
	my $densityBin = 1.0 * $tokens[1];
	my $numRadius = 1.0 * $tokens[2];
	my $radiusBin = 1.0 * $tokens[3];
	my $radiusMax = 1.0 * $tokens[4];	

	my $expectedCountAlloy = $numDensity*scalar(@elements) + $numRadius*scalar(@elements) + $numRadius*scalar(@elements)*(scalar(@elements) + 1)/2 ;
	print("Datapoints expected for alloy format: $expectedCountAlloy\n");
	my $expectedCountFS = $numDensity*scalar(@elements) + $numRadius*scalar(@elements)**2 + $numRadius*scalar(@elements)*(scalar(@elements) + 1)/2;	
	print("Datapoints expected for fs format: $expectedCountFS\n");
	
	my $count = 0; 
	
	while ($row < $lines->Count)
	{
		my $readLineFlag = 1;
		@tokens = split ' ', $lines->Item($row++);
		if (@tokens) #to ignore empty lines if there are any at the end for some reason.
		{
			for (my $i = 0; $i < scalar(@atomicNumbers); $i++)
			{
				if ($tokens[0] == $atomicNumbers[$i])
				{
					$readLineFlag = 0;
				}
			}
			if ($readLineFlag == 1) 
			{
				$count += scalar(@tokens);
			}	
		}
	} 		
	
	print ("Datapoints found: $count\n");
	if ($count == $expectedCountAlloy) 
	{
		print "Parsed alloy format.\n";
		return ALLOY;
		
	}
	elsif ($count == $expectedCountFS)
	{
	 	print "Parsed fs format.\n";
		return FS;
		
	}
	else {
		print "Error parsing format from lines.\n";
		die;
	}		
}

sub AtomicNumber
{
	my ($elementSymbol) = @_;
	my $doc = Documents->New("tmp.xsd");
	my $atom = $doc->CreateAtom($elementSymbol,Point());
	my $atomicNumber = $atom->AtomicNumber;
	$doc->Discard;
	return $atomicNumber;
} 		