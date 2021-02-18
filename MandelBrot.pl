#!perl

# Purpose: Demonstrate pure Perl running inside Materials Studio

use strict;
use warnings;

my $Cols=79; my $Lines=30;
my $MaxIter=16;
my $MinRe=-2.0; my $MaxRe=1.0;
my $MinIm=-1.0; my $MaxIm=1.0;
my @chars=(' ','.',',','-',':','/','=','H','O','A','M','%','&','$','#','@','_');

for(my $Im=$MinIm;$Im<=$MaxIm;$Im+=($MaxIm-$MinIm)/$Lines) {
    for(my $Re=$MinRe;$Re<=$MaxRe;$Re+=($MaxRe-$MinRe)/$Cols) {
        my $zr=$Re; my $zi=$Im;
        my $n;
        for($n=0;$n<$MaxIter;$n++) {
            my $a=$zr*$zr; my $b=$zi*$zi;
            if($a+$b>4.0) {
                last;
            }

            $zi=2*$zr*$zi+$Im; $zr=$a-$b+$Re;
        }
        
        print $chars[$n];
    }

    print "\n";
}
