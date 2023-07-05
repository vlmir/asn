
use strict;
use warnings;
use Carp;
my $verbose = 1;
if ( $verbose ) {
		$Carp::Verbose = 1;

##==========================================================================


sub set_metalls{

    my $acc = shift @_;

    my %metalls = ();
    $metalls{'LI'} = "lithium";  $metalls{'BE'} = "beryllium";  $metalls{'NA'} = "sodium";   $metalls{'MG'} = "magnesium";
    $metalls{'AL'} = "aluminum"; $metalls{'K'}  = "potassium";  $metalls{'CA'} = "calcium";  $metalls{'SC'} = "scandium";
    $metalls{'TI'} = "titanium"; $metalls{'V'}  = "vanadium";   $metalls{'CR'} = "chromium"; $metalls{'MN'} = "manganese";
    $metalls{'FE'} = "iron";     $metalls{'CO'} = "cobalt";     $metalls{'NI'} = "nickel";   $metalls{'CU'} = "copper";
    $metalls{'ZN'} = "zinc";     $metalls{'GA'} = "gallium";    $metalls{'RB'} = "rubidium"; $metalls{'SR'} = "strontium";
    $metalls{'Y'}  = "yttrium";  $metalls{'ZR'} = "zirconium";  $metalls{'NB'} = "niobium";  $metalls{'MO'} = "molybdenum";
    $metalls{'TC'} = "technetium"; $metalls{'RU'} = "ruthenium";$metalls{'RH'} = "rhodium";  $metalls{'PD'} = "palladium";
    $metalls{'AG'} = "silver";   $metalls{'CD'} = "cadmium";    $metalls{'IN'} = "indium";   $metalls{'SN'} = "tin";
    $metalls{'SB'} = "antimony"; $metalls{'CS'} = "caesium";    $metalls{'BA'} = "barium";   $metalls{'LA'} = "lanthanum";
    $metalls{'HF'} = "hafnium";  $metalls{'TA'} = "tantalum";   $metalls{'W'}  = "tungsten"; $metalls{'RE'} = "rhenium";
    $metalls{'OS'} = "osmium";   $metalls{'IR'} = "iridium";    $metalls{'PT'} = "platinum"; $metalls{'AU'} = "gold";
    $metalls{'HG'} = "mercury";  $metalls{'TL'} = "thallium";   $metalls{'PB'} = "lead";     $metalls{'BI'} = "bismuth";

    %{$acc} = %metalls;

}
1;
