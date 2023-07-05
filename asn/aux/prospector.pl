use strict;
use warnings;
use Carp;
my $verbose = 1;
if ( $verbose ) {
		$Carp::Verbose = 1;
	}
use Data::Dumper;

my (
$datadir, # path to data directory, no trailing '/'; R&W
@ligands, #list of ligands to screen, all upper case, optional; IN
) = @ARGV;
print "datadir:$datadir\n"; # nothing but .cif files initially !!
print "ligands: @ligands \n";
my $selected = 0;
my %prospect = ();
my %ligands = ();
my @cifs =  `ls $datadir`;
map {$ligands{$_}++} @ligands;
my $cifs = @cifs;
print "Total cif files: $cifs\n";

my $i = 0;
foreach my $cifn (@cifs){
	$i++;
	chomp $cifn;
	# printf STDERR "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%-8s %i",$cifn,$. if ($. % 100 == 0);
	my $cifpth = "$datadir/$cifn";
	open my $fhr,"< $cifpth" or croak "Can\'t open file $cifpth";

	while(<$fhr>){
		next unless (/^(ATOM   |HETATM )/);
		chomp;
		my @fs = split /\s+/;
		
		if ($ligands{$fs[2]}) {
			$prospect{$cifn}++;
			$selected++;
		}
		last if (/^#/); # block separator in cif files TODO use it for parsing iso \n
	}

	close($fhr);
	print STDERR "$i:$cifpth \n";
}
print "Read cif files: $i \n";

$i = 0;
foreach my $cifn (sort keys %prospect){
	my $cifpth = "$datadir/$cifn";
	open my $fhr,"< $cifpth" or croak "Can\'t open file $cifpth";
	$i++;
	my $ssvn = $cifn;
	substr($ssvn, -3, 3, "ssv");
	my $ssvpth = "$datadir/$ssvn";
	open my $fhw,"> $ssvpth" or croak "Failed to open file $ssvpth";
	print STDERR "$i:$ssvpth centers:$prospect{$cifn}\n";

	while(<$fhr>){
		my $sect;
		if (/^ATOM   /) {$sect = 1}
		elsif (/^HETATM /) {$sect = 2}
		else {next};
		chomp ;
		my @fs = split /\s+/;
	next if $fs[2] eq 'H';
	next if $fs[5] eq 'HOH';
		my $count = scalar @fs;
		$count =26 ?
	  print $fhw "$sect $fs[1] $fs[2] $fs[10] $fs[11] $fs[12] $fs[5] $fs[6] $fs[7]\n"
		: carp("$ssvpth: $_");
		last if (/^#/); # block separator in cif files TODO use it for parsing iso \n
	}

}
print "Saved ssv files: $i \n";
print "Selected centers: $selected \n";
