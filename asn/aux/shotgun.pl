BEGIN {
	push @INC, "/home/mironov/git/asn/t";
}
use shotgun;
use strict;
use warnings;
use Carp;
my $verbose = 1;
if ($verbose) {
$Carp::Verbose = 1;
use Data::Dumper;
};

my $start = time;
print shotgun::__date, "\n";
my (
$cifdir, # path to cif directory, no trailing '/'; R&W
$outdir,
$atmoff,
$ctrsect, # '1' || '2'
@baits, #list of baits to screen, all upper case; IN
) = @ARGV;
my $ssvdir = "$outdir/ssv";
mkdir $ssvdir unless (-e $ssvdir);
my $offdir = "$outdir/offenders";
mkdir $offdir unless (-e $offdir);
my $outpth = "$outdir/shotgun.w1";
`rm $cifdir/*.ssv $ssvdir/*.ssv`;
my %sects;
$sects{1} = "ATOM   ";
$sects{2} = "HETATM ";
croak "provide correct arg 4" unless $sects{$ctrsect} ;
my %elsy2ind = ();
my $eli = 0;
map {$elsy2ind{$_} = $eli++} sort @baits, '00';
my $pntsect = 1;
print "cifdir:$cifdir baits: @baits ball size: $atmoff\n";
my @cifs =  `ls $cifdir`;
my $cifs = @cifs;
print "Total cif files: $cifs \n";
my %preys = ();
my $i = 0; # opened cif files
my $j = 0; # found preys
my $k = 0; # cif files with preys

foreach my $cifn (@cifs){
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
			`mv $cifpth $offdir`; # excluding from subsequente re-runs
			last;
		}

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
	#rint STDERR "$i:$cifpth $sects{$ctrsect}:$j\n";
}
print "Found $j baits in $i files\n";
print shotgun::__date, "\n";
my $p = 0; # stays 0 if bounty is empty
my $dstoff = 5.0;
my %dsts = ();
my %pdbn2ind = (); # indexing in lexinumeriacal order
my $entind = 0;
my $ctrind = 0;
open my $fhw1, ">", $outpth || croak "fhw1: failed to open: $outpth";
print $fhw1 "$k $j $atmoff\n"; # file header

foreach my $pdbn (sort keys %preys){
	my %bounty = (); # now local 
	my $cifpth = "$cifdir/$pdbn.cif";
	open my $fhr,"< $cifpth" or next;
	my $ssvpth = "$ssvdir/$pdbn.ssv";
	open my $fhw,"> $ssvpth" or croak "Failed to open file $ssvpth";
	my $preys = @{$preys{$pdbn}};
	my $j = 0;
	my @ctrs = @{$preys{$pdbn}}; # array of arrays
	my @pnts = (); # array of arryas

	while(<$fhr>){
		my $ss = substr $_, 0, 7;
		next unless $ss  eq $sects{$pntsect};
		$j++;
		chomp ;
		my @fs = split /\s+/;
		next if $fs[2] eq 'H';
		next if $fs[5] eq 'HOH';
		$fs[0] = $pdbn;
			my @arr = @fs[0..2, 5..8, 10..12];
			push @{$bounty{$pdbn}}, \@arr;
		print $fhw "@arr\n"; # saving in ssv file
		last if (/^#/); # end of secton
	}
	$p++;
	close $fhw;
	#rint STDERR "$i:$ssvpth centers:$preys $sects{$pntsect}:$j\n";
	# OUTPUT
	
	$pdbn2ind{$pdbn} = $entind;
	my $ctrs = $preys{$pdbn};
	my $atms = $bounty{$pdbn};
	my $ret = shotgun::squarein($ctrs, $atms, $atmoff); # TODO fix sub calls
	my ($shots, $dsts) = @{$ret};

	foreach my $ctrid (sort {$a<=>$b} keys %{$shots}) {
		my $shot = $shots->{$ctrid};
		my $guts = @{$shot};
		carp "$pdbn:$ctrid: guts:$guts" unless $guts == $atmoff;

		foreach (@{$shot}) {
			my @gut = @{$_};
			print $fhw1 "$pdbn2ind{$pdbn} $ctrind $gut[4] @gut[0..3] $elsy2ind{$gut[-1]}\n";
		}
		$ctrind++;
	}
	$entind++;
}
close $fhw1;
carp "Unexpected diff in ctr counts: $j vs $ctrind" unless $j == $ctrind;

print "Saved ssv files: $p \n";
print "My entries: $entind My centers: $ctrind\n";
print "Range of distances to the center:\n";
print "Saved bounty in: $outpth\n";
shotgun::benchmark($start, "ALL DONE");
print shotgun::__date, "\n";

