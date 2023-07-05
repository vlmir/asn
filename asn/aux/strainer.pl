BEGIN {
	push @INC, "/home/mironov/git/asn/aux";
}
use strict;
use warnings;
use Carp;
my $verbose = 1;
if ($verbose) {
$Carp::Verbose = 1;
use Data::Dumper;
};
use strainer qw(
@allmets
$pdb
$expmaxdst
$dsts
__date
benchmark
rmsd
squarein
);
my $start = time;
print __date, "\n";
my (
$cifdir, # path to cif directory, no trailing '/'; R&W
$outdir,
$cutoff, # no filtering if <2
$ctrsect, # PDB section of the center atom, '1' || '2'
@baits, #list of baits to screen, all upper case; IN
) = @ARGV;
@baits = @allmets unless (@baits);
my $ssvdir = "$outdir/ssv";
mkdir $ssvdir unless (-e $ssvdir);
my $offdir = "$outdir/offenders";
mkdir $offdir unless (-e $offdir);
my $baitstng = join '-', @baits;
my $outpth = "$outdir/strainer-p$cutoff-$baitstng.ssv";
my %sects;
$sects{1} = "ATOM   ";
$sects{2} = "HETATM ";
croak "provide correct arg 3" unless $sects{$ctrsect} ;
my %elsy2ind = ();
my $eli = 0;
map {$elsy2ind{$_} = $eli++} sort @baits, '00';
my $pntsect = 1;
print "cifdir:$cifdir baits: @baits cloud size: $cutoff\n";
my @cifs =  `ls $cifdir`;
my $cifs = @cifs;
print "Total cif files: $cifs \n";

my %preys = ();
if ($pdb) {
### FIRST PASS ###
	my $i = 0; # opened cif files
	my $j = 0; # found preys
	my $k = 0; # cif files with preys
	foreach my $cifn (@cifs){
		## selecting entries with metals; result in %prays
		my $jj = 0; # local preys
		my $pdbn = substr $cifn, 0, 4;
		my $cifpth = "$cifdir/$pdbn.cif";
		open my $fhr,"< $cifpth" or croak "fhr: failed to open: $cifpth";
	
		while(<$fhr>){
			my $ss = substr $_, 0, 7;
			next unless $ss  eq $sects{$ctrsect};
			chomp;
			my @fs = split /\s+/;
			if (@fs != 26) {
				carp "file:$cifpth:incorrect number of fields";
				`mv $cifpth $offdir`; # excluding from subsequent re-runs
				last;
			}
			next unless $fs[25] == 1; # selecting just the 1st model
	
			if ($elsy2ind{$fs[2]}) {
				$fs[0] = $pdbn;
				my @arr = @fs[0..2, 5..8, 10..12];
				push @{$preys{$pdbn}}, \@arr;
			$j++;
			$jj++;
			}
			last if (substr $ss, 0, 1) eq '#'; # end of section
		}
		close($fhr);
		$i++;
		$k++ if $jj;
	} ## takes ~20min with the full PDB
	print "Scanned cif files: $i Found baits: $j \n";
	print __date, "\n";
	
	### SECOND PASS ###
	## extracting corresponding points and saving ssv files
	my $mycifs = 0;
	foreach my $pdbn (sort keys %preys){
		my $cifpth = "$cifdir/$pdbn.cif";
		open my $fhr,"< $cifpth" or next;
		my $ssvpth = "$ssvdir/$pdbn.ssv";
		open my $fhw,"> $ssvpth" or croak "Failed to open file $ssvpth";
		my @ctrs = @{$preys{$pdbn}}; # array of arrays
		map {print $fhw "@{$_}\n"} @ctrs; # saving in ssv file
	
		while(<$fhr>){
			my $ss = substr $_, 0, 7;
			next unless $ss  eq $sects{$pntsect}; # normally ATOM
			my @fs = split /\s+/;
			next unless $fs[25] == 1; # selecting just the 1st model
			next if $fs[2] eq 'H';
			next if $fs[5] eq 'HOH';
			$fs[0] = $pdbn;
			my @arr = @fs[0..2, 5..8, 10..12];
			print $fhw "@arr\n"; # saving in ssv file
			last if (/^#/); # end of secton
		}
		$mycifs++;
		close $fhw;
		close $fhr;
	} # ~40 min for the full PDB
} ## end of if $pdb; ~1h for the full PDB

### PARSING SSV FILES & FILTERING ###
# ~3h for the full PDB
%preys = ();
my $setpnts = 0;
my %pdbn2ind = (); # indexing in lexinumeriacal order
my $entind = 0;
my $setctrs = 0;
my $maxpnts = 0; # maximum atoms per center
my $dstmax = 0;
my $dstmin = $expmaxdst;

my @ssvs =  `ls $ssvdir`;
open my $fhw1, ">", $outpth || croak "fhw1: failed to open: $outpth";
foreach my $ssvn (@ssvs){
	my %bounty = (); # now local 
	my $pdbn = substr $ssvn, 0, 4;
	my $ssvpth = "$ssvdir/$pdbn.ssv";
	open my $fhr,"< $ssvpth" or croak "fhr: failed to open: $ssvpth";

	while(<$fhr>){
		my @fs = split /\s+/;
		if (@fs != 10) {
			croak "pdbn:$pdbn fhr: wrong number of fields"
		}
		if ($elsy2ind{$fs[2]}) {
			push @{$preys{$pdbn}}, \@fs;
		} else {
			push @{$bounty{$pdbn}}, \@fs;
		}
	}
	close($fhr);
	
	$pdbn2ind{$pdbn} = $entind; # for mapping pdbn <-> entind
	my $ctrs = $preys{$pdbn}; # centers in the entry
	my $atms = $bounty{$pdbn}; # ALL atoms in the entry
	my $shots = squarein($ctrs, $atms, $cutoff);

	foreach my $ctrid (sort {$a<=>$b} keys %{$shots}) {
		my $shot = $shots->{$ctrid};
		my $count = @{$shot}; # number of atoms in the center
		## the line below makes sense only for a fixed number of atoms per center
		carp "warning: $pdbn:$ctrid: atoms:$count != cutoff:$cutoff" unless $count == $cutoff;
		$maxpnts = $count if $count > $maxpnts;
		foreach (@{$shot}) {
			my @gut = @{$_};
			my $dst = $gut[3];
			print $fhw1 
			"$setctrs @gut[0..3] $elsy2ind{$gut[-1]} $gut[4] $pdbn2ind{$pdbn}\n";
			$dstmax = $dst if $dst > $dstmax;
			$dstmin = $dst if ($dst < $dstmin && $dst > 0);
			$setpnts++;
		}
		$setctrs++;
	}
	$entind++;
} # end of for each ssvs
close $fhw1;

print "Distribution maxdst/ctr: count \n";
map {print $_-1, "..$_ A: $dsts->{$_} \n"} sort {$a<=>$b} keys %{$dsts};
print "MyCenters: $setctrs MaxPoints: $maxpnts MyPoints: $setpnts \n";
print "Range of distances to the center: $dstmin..$dstmax \n";
print "Saved bounty in: $outpth\n";
benchmark($start, "ALL DONE");
print __date, "\n";

