## Tetrahedral.pm

package Tetrahedral;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 4
			  );

our $expectedAngle = 109.5;

our $invCorrM = [
	[0.722, -0.111, -0.111, -0.111, -0.111, -0.278],
	[-0.111, 0.722, -0.111, -0.111, -0.278, -0.111],
        [-0.111, -0.111, 0.722, -0.278, -0.111, -0.111],
        [-0.111, -0.111, -0.278, 0.722, -0.111, -0.111],
        [-0.111, -0.278, -0.111, -0.111, 0.722, -0.111],
        [-0.278, -0.111, -0.111, -0.111, -0.111, 0.722]
];

sub new
  {
  my $class = shift @_;
  my $self = Coordination->new(@defaultDataMembers, @_);

  bless $self, ref $class || $class;

  return $self;
  }


sub orderedCombinations
  {
  my $self = shift @_;

  return ([@_]);
  }

sub angleTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $combo = shift @_;
  my $angleStats = shift @_;
  my $distChi = shift @_;
  my $leaveOut = (@_)? shift @_: 0;
 
  my ($mean109, $varianceOrN);

  if ($type eq "dev")
    {
    $varianceOrN = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle) ** 2) / $varianceOrN ;} ($self->calcAllAngles109($combo)) ;

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean109 = ($$angleStats{"tetrahedral"}{"109"}{"mean"})? ($$angleStats{"tetrahedral"}{"109"}{"mean"}) : $expectedAngle;
      #$varianceOrN = ($$angleStats{"tetrahedral"}{"109"}{"variance"})? ($$angleStats{"tetrahedral"}{"109"}{"variance"}) : $$angleStats{"variance"}; ## major variance
      $varianceOrN = $$angleStats{"variance"}; ## overall variance
      }
    elsif ($type eq "ownStats")
      {
      $mean109 = $$angleStats{"tetrahedral"}{"109"}{"mean"};
      $varianceOrN = $$angleStats{"tetrahedral"}{"109"}{"variance"};
      }

    #my @means = ($expectedAngle, $expectedAngle, $expectedAngle, $expectedAngle, $expectedAngle, $expectedAngle); ## ideal mean
    my @means = ($mean109, $mean109, $mean109, $mean109, $mean109, $mean109); ## major mean
    my @angles = $self->calcAllAngles109($combo);
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    ## set the smallest angle's diff to zero
    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std = 1/sqrt($varianceOrN);
    my $invStds = [$std, $std, $std, $std, $std, $std];

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_])**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob0 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.1)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob1 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.2)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob2 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.3)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob3 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.4)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob4 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.5)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob5 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.6)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob6 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.7)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob7 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.8)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob8 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.9)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob9 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 2)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob10 = &Statistics::Distributions::chisqrprob(10, $chiStat);

    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);
    my $prob = &Statistics::Distributions::chisqrprob(9, $chiStat+ $distChi);
print $self->{shellObj}->metalID(), ", ", ref $self, ", $prob0, $prob1, $prob2, $prob3, $prob4, $prob5, $prob6, $prob7, $prob8, $prob9, $prob10, $prob\n" if $prob>0.5;

    return $chiStat;
    }
  }


sub angleList
  {
  my $self = shift @_;

  return 0 if (! exists $self->{bestCombo});

  my @anglelist;
  push (@anglelist, $self->calcAllAngles109());

  my $pair1 = [sort {$a <=> $b} ($anglelist[0], $anglelist[5])];
  my $pair2 = [sort {$a <=> $b} ($anglelist[1], $anglelist[4])];
  my $pair3 = [sort {$a <=> $b} ($anglelist[2], $anglelist[3])];

  my @pairlist = sort { $$a[0] <=> $$b[0]; } ($pair1, $pair2, $pair3);
  @anglelist = ( @{$pairlist[0]}, @{$pairlist[1]}, @{$pairlist[2]} );

  return @anglelist;
  }

## ## Order angles as largest, middles, and opposite, with the ligands ordered as
#	first is the ligand that composing the largest, and is also the first ligand defining the largest angle plane where the perpendicular vector and one of the rest two is the closest.
#	second is the other of the largest-angle ligand,
#	third is the one closest to the perpendicular vector
#	forth is the last ligand.

sub orderedAnglesPerpVec
  {
  my $self = shift @_;
  return 0 if (! exists $self->{bestCombo});
  my $combo = $self->{bestCombo}->{ligands};

  my @anglelist;
  push (@anglelist, $self->calcAllAngles109());
  @anglelist = sort {$b <=> $a} (@anglelist);

  my $center = $self->{shellObj}->{center};
  my ($aa, $bb);
  FORLOOP: for(my $x = 0; $x < $#$combo; $x++)
    {
    for(my $y = $x+1; $y < @$combo; $y++)
      {
      my $angle = $center->angle($$combo[$x], $$combo[$y]);

      if ($angle == $anglelist[0])
        {
        $aa = $$combo[$x];
        $bb = $$combo[$y];
	last FORLOOP;
        }
      }
    }

  my @abv = $center->perpendicularVec($aa, $bb);
  my @bav = $center->perpendicularVec($bb, $aa);

  my $abp_x = $abv[0] + $center->{x};
  my $abp_y = $abv[1] + $center->{y};
  my $abp_z = $abv[2] + $center->{z};
  my $bap_x = $bav[0] + $center->{x};
  my $bap_y = $bav[1] + $center->{y};
  my $bap_z = $bav[2] + $center->{z};

  my $abp = Atom->new("x" => $abp_x, "y" => $abp_y, "z" => $abp_z);
  my $bap = Atom->new("x" => $bap_x, "y" => $bap_y, "z" => $bap_z);

  my @cd = grep { $_ ne $aa && $_ ne $bb;} (@$combo);
  my $angleAbpcd0 = $center->angle($abp, $cd[0]);
  my $angleAbpcd1 = $center->angle($abp, $cd[1]);
  my $angle_bapcd0 = $center->angle($bap, $cd[0]);
  my $angle_bapcd1 = $center->angle($bap, $cd[1]);

  my $minAngle = (sort {$a <=> $b} ($angleAbpcd0, $angleAbpcd1, $angle_bapcd0, $angle_bapcd1))[0];
  my $newCombo = [];
  if ($angleAbpcd0 == $minAngle) 
    { @$newCombo = ($aa, $bb, $cd[0], $cd[1]);  }
  elsif ($angleAbpcd1 == $minAngle)
    { @$newCombo = ($aa, $bb, $cd[1], $cd[0]);  }
  elsif ($angle_bapcd0 == $minAngle)
    { @$newCombo = ($bb, $aa, $cd[0], $cd[1]);  }
  elsif ($angle_bapcd1 == $minAngle)
    { @$newCombo = ($bb, $aa, $cd[1], $cd[0]);  }

  $self->{bestCombo}->{ligands} = $newCombo;
  return $self->calcAllAngles109();
  }


## Order angles as largest, middles, and opposite, with the ligands ordered as 
#	the ligand sharing between largest and largest-middle first, 
#	the other ligand composing the largest second, 
#	the other ligand composing the largest-middle third,
#	and the fourth. 

sub orderedAnglesLargeMid
  {
  my $self = shift @_;
  return 0 if (! exists $self->{bestCombo});
  my $combo = $self->{bestCombo}->{ligands};

  my @anglelist;
  push (@anglelist, $self->calcAllAngles109());
  @anglelist = sort {$b <=> $a} (@anglelist);
  
  my $center = $self->{shellObj}->{center};
  my (%first, %second, %third);
  for(my $x = 0; $x < $#$combo; $x++)
    {
    for(my $y = $x+1; $y < @$combo; $y++)
      {
      my $angle = $center->angle($$combo[$x], $$combo[$y]);

      if ($angle == $anglelist[0]) 
	{ 
	$first{$x} = 1;
	$first{$y} = 1;
	}
      elsif ($angle == $anglelist[1])
        {
	$second{$x} = 1;
        $second{$y} = 1;
        }
      elsif ($angle == $anglelist[2])
        {
	$third{$x} = 1;
	$third{$y} = 1;
	}
      }
    }

  unless ( grep {$second{$_} ;} (keys %first) ) {%second = %third};
 
  my @ind;
  push @ind, grep {$second{$_} == 1;} (keys %first);
  push @ind, grep {$second{$_} != 1;} (keys %first);
  push @ind, grep {$first{$_} != 1;} (keys %second);
  push @ind, grep {! ($first{$_} == 1 || $second{$_} == 1);} (0..3);

  my $newCombo = [];
  map {$$newCombo[$_] = $$combo[$ind[$_]];} (0..3);
  $self->{bestCombo}->{ligands} = $newCombo;

  return $self->calcAllAngles109();
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"tetrahedral"}{"109"}}, $_->calcAllAngles109());
  
  return 0;
  }

sub calcAllAngles109
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};
  
  my $center = $self->{shellObj}->{center};

  my @angles;
  for(my $x = 0; $x < $#$combo; $x++)
    {
    for(my $y = $x+1; $y < @$combo; $y++)
      {
      push (@angles, $center->angle($$combo[$x], $$combo[$y]));
      }
    }

  return @angles;
  }




