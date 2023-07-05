package shotgun;
use strict;
use warnings;
use Carp;
my $verbose = 1;
if ( $verbose ) {
	$Carp::Verbose = 1;
}
use Data::Dumper;
my $expmaxdst = 500;

sub __date {
        my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
         # e.g. 2008:05:11 12:52
        my $result = sprintf "%4d:%02d:%02d %02d:%02d", $year+1900,$mon+1,$mday,$hour,$min;
}

sub benchmark {
        my (
        $start_time,
        $message,
        $date, # boolean, optional
        ) = @_;
        my $elapsed_time = (time - $start_time)/60;
        $message = "$message in ";
        $date ?
        print $message, sprintf("%.2f", $elapsed_time), " min ", __date(), "\n" :
        print $message, sprintf("%.2f", $elapsed_time), " min\n";
}

sub rmsd {
	my (
		$pnt1,
		$pnt2
	) = @_;
	my $ss = 0;

	for (my $i = 0; $i < 3; $i++) {
		$ss += ($pnt1->[$i] - $pnt2->[$i])*($pnt1->[$i] - $pnt2->[$i]);
	}
	return sqrt($ss);
}

sub squarein {
	my (
		$arr1,
		$arr2,
		$cutoff
	) = @_;
	my $dsts; # {}[];
	my $shots; # {}[][]
	foreach my $ctr (@{$arr1}) {
		my $ctratms = 0;
		my $ctrid = $ctr->[1];
		my $ctrelsy = $ctr->[2];
		my @ctrxyz = @{$ctr}[7..9];
		my $shot; # [][]
		my %ind2dst = ();
		my @myself = ( @ctrxyz, 0.0, $ctrid, $ctrelsy);
		push @{$shot}, \@myself;
		my $atms = 0; # the total number of atoms in the end, not index!

		# computing all distances
		foreach my $atm (@{$arr2}) {
			my @atmxyz = @{$atm}[7..9];
			my $dst = shotgun::rmsd(\@ctrxyz, \@atmxyz);
			$ind2dst{$atms} = $dst;
			$atms++;
		}
		# selecting atoms closest to the center
		my @sdsts = sort {$a<=>$b} values %ind2dst;
		$dsts->{max} = $sdsts[-1];
		$dsts->{min} = $sdsts[0];
		my $dstoff = $sdsts[$cutoff-2]; # dstoff excludes ctr atom
		$dsts->{off} = $dstoff;
		my %myi2d;
		foreach my $ind (keys %ind2dst) {
			my $dst = $ind2dst{$ind};
			if ($dst <= $dstoff) {
				$myi2d{$ind} = $dst;
			}
		}

		foreach my $ind (sort {$a<=>$b} keys %myi2d) {
			my $atm = $arr2->[$ind];
			my @buddies = (@{$atm}[7..9], $ind2dst{$ind}, $atm->[1], '00');
			push @{$shot}, \@buddies;
		}
		$shots->{$ctrid} = $shot;
	}
	my @ret = ($shots, $dsts);
	return \@ret;
}
###############################################################################
sub roundup {
	my (
		$arr1,
		$arr2,
$cutoff
	) = @_;
my $newshots; # {}[][]
	my $mindst = $cutoff;
	my $maxdst = 0.0;
		foreach my $arr1 (@{$arr1}) {
		my $atms = 0;
		my $ctrid = $arr1->[1];
		my @ctrxyz = @{$arr1}[7..9];
		push @ctrxyz, 0.0;
		my $preys = $newshots->{$ctrid}; # [][]
		$preys->[0] = \@ctrxyz;

		foreach my $arr2 (@{$arr2}) {
		my @atmxyz = @{$arr2}[7..9];
		my $dst = shotgun::rmsd(\@ctrxyz, \@atmxyz);
		next unless ($dst <= $cutoff);
		push @atmxyz, $dst;
		push @{$preys}, \@atmxyz;
		$atms++;
		$mindst = $dst if ($dst < $mindst);
		$maxdst = $dst if ($dst > $maxdst);
		}
		my $xyzsz = @{$preys};
		print STDERR "ctrid:$ctrid atms:$atms mindst:$mindst maxdst:$maxdst\n";
	}
	return $newshots;
}
1;
