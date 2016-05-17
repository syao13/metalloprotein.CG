## Octahedral.pm 

package Octahedral;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 6
			  );

our $expectedAngle90 = 90;
our $expectedAngle180 = 180;

our $invCorrM = [
	[0.625,	-0.125,	-0.125,	-0.375,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000],
	[-0.125, 0.625,	-0.375,	-0.125,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000],
	[-0.125,-0.375,	0.625,	-0.125,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000],
	[-0.375,-0.125,	-0.125,	0.625,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000],
	[0.000,	0.000,	0.000,	0.000,	0.625,	-0.125,	-0.125,	-0.375,	0.000,	0.000,	0.000,	0.000],
	[0.000,	0.000,	0.000,	0.000,	-0.125,	0.625,	-0.375,	-0.125,	0.000,	0.000,	0.000,	0.000],
	[0.000,	0.000,	0.000,	0.000,	-0.125,	-0.375,	0.625,	-0.125,	0.000,	0.000,	0.000,	0.000],
	[0.000,	0.000,	0.000,	0.000,	-0.375,	-0.125,	-0.125,	0.625,	0.000,	0.000,	0.000,	0.000],
	[0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.625,	-0.125,	-0.125,	-0.375],
	[0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	-0.125,	0.625,	-0.375,	-0.125],
	[0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	-0.125,	-0.375,	0.625,	-0.125],
	[0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	-0.375,	-0.125,	-0.125,	0.625]
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
  my $combo = shift @_;

  my $orderedCombos;
  foreach my $oppositeIndex (1..$#$combo)
    {
    my @restIndeces = grep { $oppositeIndex != $_; } (1..$#$combo);
    push (@$orderedCombos, [$$combo[0], $$combo[$oppositeIndex], map { $$combo[$restIndeces[$_]]; } (0,1,2,3)] );
    push (@$orderedCombos, [$$combo[0], $$combo[$oppositeIndex], map { $$combo[$restIndeces[$_]]; } (0,2,1,3)] ); 
    push (@$orderedCombos, [$$combo[0], $$combo[$oppositeIndex], map { $$combo[$restIndeces[$_]]; } (0,3,1,2)] );
    }

  return $orderedCombos;
  }

sub angleTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $combo = shift @_;
  my $angleStats = shift @_;
  my $distChi = shift @_;
  my $leaveOut = (@_)? shift @_: 0;

  my ($mean90, $varianceOrN90, $varianceOrN180);

  if ($type eq "dev")
    {
    $varianceOrN90 = $self->numAngles() ;
    $varianceOrN180 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle90) ** 2) / $varianceOrN90 ;} ($self->calcAllAngles90($combo)) ; #ideal mean, overal variance
    map { $angleTestStat += (($_ - $expectedAngle180) ** 2) / $varianceOrN180 ;} ($self->calcAllAngles180($combo)) ;

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean90 = ($$angleStats{"octahedral"}{"90"}{"mean"})? ($$angleStats{"octahedral"}{"90"}{"mean"}) : $expectedAngle90;
      #$varianceOrN90 = ($$angleStats{"octahedral"}{"90"}{"variance"})? ($$angleStats{"octahedral"}{"90"}{"variance"}) : $$angleStats{"variance"}; ## major coordination variance
      $varianceOrN90 = $$angleStats{"variance"}; ## overall variance
      }
    else
      {
      $varianceOrN90 =  $$angleStats{"octahedral"}{"90"}{"variance"};
      }

    #my @means = ($expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90); ## ideal mean
    my @means = ($mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90, $mean90); ## major coordinaiton mean
    my @angles = ($self->calcAllAngles90($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std90 = 1/sqrt($varianceOrN90);
    my $invStds = [$std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90, $std90];

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_])**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob0 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.1)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob1 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.2)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob2 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.3)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob3 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.4)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob4 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.5)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob5 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.6)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob6 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.7)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob7 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.8)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob8 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 1.9)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob9 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = 0;
    map {$chiStat += ($$diff[$_] * $$invStds[$_] / 2)**2;} (0..(@$diff-1));
    $chiStat += $distChi;
    my $prob10 = &Statistics::Distributions::chisqrprob(21, $chiStat);

    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);
    my $prob = &Statistics::Distributions::chisqrprob(15, $chiStat + $distChi);
print $self->{shellObj}->metalID(), ", ", ref $self, ", $prob0, $prob1, $prob2, $prob3, $prob4, $prob5, $prob6, $prob7, $prob8, $prob9, $prob10, $prob\n" if $prob>0.5;

    #print "mean 90, $expect90\nangles: ";
    #print map {"$_, "; } (@angles);
    #print "\nstd90, $std90\n";
    #print "covariance: ", $chiStat;#, "; previous: $angleTestStat\n";

    return $chiStat;
    }
  }


sub angleList
  {
  my $self = shift @_;

  return 0 if (! exists $self->{bestCombo});

  my @anglelist;
  push (@anglelist, $self->calcAllAngles90());
  push (@anglelist, $self->calcAllAngles180());

  return @anglelist;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"octahedral"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"octahedral"}{"180"}}, $_->calcAllAngles180());

  return 0;
  }


sub calcAllAngles90
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my $center = $self->{shellObj}->{center};

  my @angles;
  my $pair1 = [$$combo[0], $$combo[1]];
  my $pair2 = [$$combo[2], $$combo[3]];
  my $pair3 = [$$combo[4], $$combo[5]];

  for(my $x = 0; $x < @$pair1; $x++)
    {
    for(my $y = 0; $y < @$pair2; $y++)
      {
      push (@angles, $center->angle($$pair1[$x], $$pair2[$y])) ;
      }
    }

  for(my $x = 0; $x < @$pair1; $x++)
    {
    for(my $y = 0; $y < @$pair3; $y++)
      {
      push (@angles, $center->angle($$pair1[$x], $$pair3[$y]));
      }
    }

  for(my $x = 0; $x < @$pair2; $x++)
    {
    for(my $y = 0; $y < @$pair3; $y++)
      {
      push (@angles, $center->angle($$pair2[$x], $$pair3[$y])) ;
      }
    }

  return @angles;
  }

sub calcAllAngles180
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my $center = $self->{shellObj}->{center};

  my @angles;
  my $pair1 = [$$combo[0], $$combo[1]];
  my $pair2 = [$$combo[2], $$combo[3]];
  my $pair3 = [$$combo[4], $$combo[5]];

  foreach my $pair ($pair1, $pair2, $pair3)
    {
    for(my $x = 0; $x < $#$pair; $x++)
      {
      for(my $y = $x+1; $y < @$pair; $y++)
        {
        push (@angles, $center->angle($$pair[$x], $$pair[$y])) ;
        }
      }
    }

  return @angles;
  }













