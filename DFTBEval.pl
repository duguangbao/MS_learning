
#!perl

use strict;
use MaterialsScript qw(:all);

#############################################################
#####################  DESCRIPTION  #########################
#############################################################
# This script will run DFTB+ and DMol3 geometry optimization  
# calculations on all .xsd documents in the launch folder and  
# compare any distance measurements and angle measurements
# defined in the documents. It will also calculate and compare 
# atomization energy for any non periodic structures.
# Bond lengths are given in Angstrom and Angles in degree.
# The atomization difference is given in kcal/mol.
# At the end of the output file statistics for the distances
# and angles will be reported.
#############################################################
######################  SETTINGS  ###########################
#############################################################


# Input name of .skflib file.
# In the case of a default library input the name of the library.
my $skflibFile = "CHNO";

# Input settings for DMol3 and DFTB.
# DMol3 functional settings should be the same as was used for the parameterization
my $DMol3Settings = Settings(
                    Quality => "Fine",
                    TheoryLevel => "GGA",
                    NonLocalFunctional => "PBE",
                    BasisFile => "4.4",
                    UseSmearing => "Yes",
                    OptimizeCell => "Yes"
                     );
my $DFTBSettings = Settings(
                    Quality => "Fine",
                    UseSmearing => "Yes",
                    SKFLibrary => $skflibFile,
                    OptimizeCell => "Yes",
                    WriteLevel => "Silent"
                     );       		     
                   

#############################################################
############## NO CHANGES BEYOND THIS LINE###################
#############################################################


Documents->SaveAllAtEnd = "Yes";

my $output = Documents->New("output.txt");
my %doclist = %Documents;
my %distNumber = ("All",0);
my %distAbsDiff = ("All",0.0);
my %anglNumber = ("All",0);
my %anglAbsDiff = ("All",0.0);
my %atomEnergyDmol3 = ();
my %atomEnergyDFTB = ();
my $AEnergyDMol3;
my $AEnergyDFTB;

#Loop over all documents
foreach my $file (sort{$doclist{$a}->Name cmp $doclist{$b}->Name} sort keys %doclist)
{
    my $dist1;
    my $dist2;
    my $angl1;
    my $angl2;
    my $doc2;
    my $doc = $Documents{$file};
    
    #Only evaluate atomistic documents
    if( $doc->Type eq "3DAtomistic" )
    {
        #Number the atoms for identification
        my $number = 0;
        my $atoms;
        if (($doc->SymmetrySystems->Count == 0) || ($doc->SymmetrySystem->SymmetryDefinition->Periodicity == 0)){
            $atoms = $doc->Atoms;
        }
        else {
            $atoms = $doc->UnitCell->Atoms;
        }

        foreach my $atom (@$atoms)
        {
            $number++;
            my $name = $atom->Name;
            $name =~s/([A-Z])\d*/$1/;
            $name .=$number;
            $atom->Name = $name;
        }

        #Output document name, output details will folllow
        $output->Append($doc->Name."\n");
	my $underline = "-" x length($doc->Name);
	$output->Append($underline."\n");

        #DFTB+ and DMol3 geometry optimizations
        my $moduleName;
        my $err = "False";
        eval{
            #Any failure is caught and reported
          
            #DFTB+ geometry optimizations
            $moduleName = "DFTB+";
            my $DFTBData = Modules->DFTB->GeometryOptimization->Run($doc, $DFTBSettings );
            $DFTBData->DFTBReport->Delete;

            #Calculate atomization energy if  structure is non periodic
            my $DFTBModule = Modules->DFTB;
            if (($doc->SymmetrySystems->Count == 0) || ($doc->SymmetrySystem->SymmetryDefinition->Periodicity == 0))
            {
                $AEnergyDFTB = calculateAtomizationEnergy($DFTBModule, $DFTBSettings, $doc, \%atomEnergyDFTB);
            }
            $dist2 = $doc->UnitCell->Distances;
            $angl2 = $doc->UnitCell->Angles;
          
            #Create document copy for DMol3
            $doc2 = Documents->New($doc->Name.".xsd");
            $doc2->CopyFrom($doc);
          
            #Dmol3 geometry optimizations
            $moduleName = "DMol3";
            my $DMol3Data = Modules->DMol3->GeometryOptimization->Run($doc2, $DMol3Settings );
            $DMol3Data->Report->Delete;

            #Calculate atomization energy if  structure is non periodic
            my $DMOl3Module = Modules->DMol3;
            if (($doc->SymmetrySystems->Count == 0) || ($doc->SymmetrySystem->SymmetryDefinition->Periodicity == 0))
            {
                $AEnergyDMol3 = calculateAtomizationEnergy($DMOl3Module, $DMol3Settings, $doc, \%atomEnergyDmol3, $DMol3Data);
            }

            $dist1 = $doc2->UnitCell->Distances;
            $angl1 = $doc2->UnitCell->Angles;
        
            #Return true, no need for error handling
            1;
        } or do{
            #Report a failure
            $err = "True";
            $output->Append("Calculation failed using $moduleName!\n$@ \n");
            $output->Append("==============================================\n\n");
            $output->Save;
        };
    
        #Compare geometries
        if ($err eq "False")
        {
            #Ensure that we have any distances to report
            if ($dist1->Count > 0 )
            {
            #Output distances from DMol3 and DFTB+
            &outputDistances($output, $dist1, $dist2);
            &collectDistanceStatistics($dist1, $dist2, \%distNumber, \%distAbsDiff);
            }
        
            #Ensure that we have any angles to report
            if ($angl1->Count > 0)
            {
                #Output angles from DMol3 and DFTB+
                &outputAngles($output, $angl1, $angl2);
                &collectAngleStatistics($angl1, $angl2, \%anglNumber, \%anglAbsDiff);
            }
            
            if (($doc->SymmetrySystems->Count == 0) || ($doc->SymmetrySystem->SymmetryDefinition->Periodicity == 0))
            {
                my $AtomizationDiff = $AEnergyDMol3 - $AEnergyDFTB;
                my $AtomizationDiff = sprintf "%.5f", $AtomizationDiff;
                $output->Append("\nAtomization Diff = $AtomizationDiff\n");
            }
            $output->Append("==============================================\n");
            $output->Append("\n");
            $output->Save;
        }
    }
    $doc->Close;
}

#Output distance error statisitics
outputStatistics($output, "Bond Error Statistics:",\%distNumber, \%distAbsDiff);
#Output angle error statisitics
outputStatistics($output, "Angle Error Statistics:",\%anglNumber, \%anglAbsDiff);

######################################################
################  END OF MAIN  #######################
######################################################

#Output distances and distances differences for DMol3 and DFTB+
sub outputDistances
{
    my ($output, $dist1, $dist2) = @_;
    #Output distances from DMol3 and DFTB+
    $output->Append("DMol3 ");
    foreach my $d (@$dist1)
    {
        my $val = sprintf "%-8.5f", $d->Distance;
        $output->Append($d->Name." = ".$val." ");
    }
    $output->Append("\n");

    $output->Append("DFTB+ ");
    foreach my $d (@$dist2)
    {
        my $val = sprintf "%-8.5f", $d->Distance;
        $output->Append($d->Name." = ".$val." ");
    }
    $output->Append("\n");

    $output->Append("Diff  ");
    for (my $index = 0; $index < $dist2->Count; $index++)
    {
        my $diff = $dist2->Item($index)->Distance - $dist1->Item($index)->Distance;
        my $diffval = sprintf "%-8.5f", $diff;
        #Output angle differences between DMol3 and DFTB+
        $output->Append($dist1->Item($index)->Name." = $diffval ");
    }
    $output->Append("\n");
}

#Output angles and angel differences for DMol3 and DFTB+
sub outputAngles
{
    my ($output, $angl1, $angl2) = @_;
    $output->Append("\n");
    #Output angles from DMol3 and DFTB+
    $output->Append("DMol3 ");
    foreach my $a (@$angl1)
    {
        my $val = sprintf "%-10.5f", $a->Angle;
        $output->Append($a->Name." = ".$val." ");
    }
    $output->Append("\n");

    $output->Append("DFTB+ ");
    foreach my $a (@$angl2)
    {
        my $val = sprintf "%-10.5f", $a->Angle;
        $output->Append($a->Name." = ".$val." ");
    }
    $output->Append("\n");

    $output->Append("Diff  ");
    for (my $index = 0; $index < $angl2->Count; $index++)
    {
        my $diff = $angl2->Item($index)->Angle-$angl1->Item($index)->Angle;
        my $diffval = sprintf "%-10.5f", $diff;
        #Output angle differences between DMol3 and DFTB+
        $output->Append($angl1->Item($index)->Name." = $diffval ");
    }
    $output->Append("\n");
}

#Collect data for distance statistics
sub collectDistanceStatistics
{
    my ($dist1, $dist2, $distNumber, $distAbsDiff) = @_;
  
    for (my $index = 0; $index < $dist2->Count; $index++)
    {
        #Save data for statistics
        my $d1 = $dist1->Item($index);
        my $d2 = $dist2->Item($index);
        my $diff = $d2->Distance-$d1->Distance;
        my $tag = $d1->Name;
        
        #Ensure that data is saved under a singel tag either X-Y or Y-X
        $tag =~s/([A-Z])\d+-([A-Z])\d+/$1-$2/;
        my $tag2 = $2."-".$1;
        if($distNumber->{$tag} || $distNumber->{$tag2})
        {
            $tag = $tag2 if($distNumber->{$tag2});
            $distNumber->{$tag}++;
            $distAbsDiff->{$tag} += abs($diff);
        }
        else
        {
            $distNumber->{$tag} = 1;
            $distAbsDiff->{$tag} = abs($diff);
        }
        #Total error statistics
        $distNumber->{"All"}++;
        $distAbsDiff->{"All"} += abs($diff);
    }
}
#Collect data for angle statistics
sub collectAngleStatistics
{
    my ($angl1, $angl2, $anglNumber, $anglAbsDiff) = @_;
  
    for (my $index = 0; $index < $angl2->Count; $index++)
    {
        #Save data for statistics
        my $d1 = $angl1->Item($index);
        my $d2 = $angl2->Item($index);
        my $diff = $d2->Angle-$d1->Angle;
        my $tag = $d1->Name;
        
        #Ensure that data is saved under a singel tag either X-Y-Z or Z-Y-X
        $tag =~s/([A-Z])\d+-([A-Z])\d+-([A-Z])\d+/$1$2$3/;
        my $tag2 = $3.$2.$1;
        if($anglNumber->{$tag} || $anglNumber->{$tag2})
        {
            $tag = $tag2 if($anglNumber->{$tag2});
            $anglNumber->{$tag}++;
            $anglAbsDiff->{$tag} += abs($diff);
        }
        else
        {
            $anglNumber->{$tag} = 1;
            $anglAbsDiff->{$tag} = abs($diff);
        }
        #Total error statistics
        $anglNumber->{"All"}++;
        $anglAbsDiff->{"All"} += abs($diff);
    }
}

#Output statistics first atom specific and then total average
sub outputStatistics
{
    my ($output, $title, $number, $absDiff) = @_;

    $output->Append("\n");
    $output->Append("$title\n");
    foreach my $tag (sort keys %$number)
    {
        if($tag ne "All")
        {
            my $mean = $absDiff->{$tag}/$number->{$tag};
            $mean = sprintf "%.5e", $mean;
            $output->Append("$tag = $mean\n");
        }
    }
    $output->Append("=================\n");
    my $totalmean;
    if($anglNumber{"All"} > 0){
        $totalmean = $absDiff->{"All"}/$number->{"All"};
        $totalmean = sprintf "%.5e", $totalmean;
        $output->Append("Total Average = $totalmean\n");
    }
}

###########################################
#Atomization energy utility below this line
###########################################

#Calculate atomization energy for a structure
#defined in $doc. $atomEnergy is passed back
#to avoid calculating the energy for an atom
#more than once.
#The module and settings are passed in so that
#the same code can be used for DMol3 and DFTB+
sub calculateAtomizationEnergy
{
    my ($module, $settings, $doc, $atomEnergy, $result) = @_;
    my %atomList=();
    makeAtomList($doc,\%atomList);


    updateAtomEnergyList($module, $settings,\%atomList,$atomEnergy);
    my $totalEnergyAtoms = 0.0;
    foreach my $key (sort keys %atomList) {
        $totalEnergyAtoms += $atomEnergy->{$key}*$atomList{$key};
    }
    my $totalEnergyMolecule;
  
    if($result)
    {
        $totalEnergyMolecule = $result->TotalEnergy;
    } 
    else
    {
        $totalEnergyMolecule = $doc->PotentialEnergy;
    }
    
    return $totalEnergyMolecule - $totalEnergyAtoms;	
}

#Calculate atomic energies if data does not already exist
sub updateAtomEnergyList
{
    my ($module, $settings, $atomList, $atomEnergy) = @_;
    foreach my $key (sort keys %{$atomList}) {
        if(! exists $atomEnergy->{$key})
        {
            $atomEnergy->{$key} = makeSingleAtomCalulation($module,$settings,$key);
    
        }
    }
}

#Create list of unique elements and how many they are
sub makeAtomList
{
    my ($doc, $atomList) = @_;

    my $atoms = $doc->Atoms;

    foreach my $atom (@$atoms) {
        if(exists $atomList->{$atom->ElementSymbol})
        {
            $atomList->{$atom->ElementSymbol} += 1;
        }
        else
        {
            $atomList->{$atom->ElementSymbol} = 1;
        }
    }
}

#Calculate total energy for a single atom
#The module and settings are passed in so that
#the same code can be used for DMol3 and DFTB+
sub makeSingleAtomCalulation
{
    my ($module,$settings,$element) = @_;

    my $newdoc;
    $newdoc= Documents->New($element.".xsd");
    $newdoc->CreateAtom( $element, Point(X => 0.0, Y => 0.0, Z => 0.0) );

    my $result = $module->Energy->Run($newdoc,$settings);
    my $TotalEnergy;
  
    $TotalEnergy = $newdoc->PotentialEnergy;
    if(! $TotalEnergy)
    {
        $TotalEnergy = $result->TotalEnergy;
    }
    eval{$result->Report->Delete;};
    eval{$result->DFTBReport->Delete;};
    $newdoc->Delete;
    return $TotalEnergy;
}
