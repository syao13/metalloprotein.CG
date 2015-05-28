## SquarePyramidalV.pm 

package SquarePyramidalV;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 4
			  );

our $expectedAngle90a = 90;
our $expectedAngle90p = 90;
our $expectedAngle180 = 180;

our $invCorrM = [
	[1,	0.0000000,	0.0000000,	0.0000000,	0.0000000],
	[0,	1.3333333,	0.6666667,	0.0000000,	0.0000000],
	[0,	0.6666667,	1.3333333,	0.0000000,	0.0000000],
	[0,	0.0000000,	0.0000000,	1.3333333,	0.6666667],
	[0,	0.0000000,	0.0000000,	0.6666667,	1.3333333]
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
  push (@$orderedCombos, [$$combo[0], $$combo[1], $$combo[2], $$combo[3]]);
  push (@$orderedCombos, [$$combo[0], $$combo[2], $$combo[1], $$combo[3]]);
  push (@$orderedCombos, [$$combo[0], $$combo[3], $$combo[1], $$combo[2]]);
  push (@$orderedCombos, [$$combo[2], $$combo[3], $$combo[0], $$combo[1]]);
  push (@$orderedCombos, [$$combo[1], $$combo[3], $$combo[0], $$combo[2]]);
  push (@$orderedCombos, [$$combo[1], $$combo[2], $$combo[0], $$combo[3]]);

  return $orderedCombos;
  }

sub angleTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $combo = shift @_;
  my $angleStats = shift @_;
  my $leaveOut = (@_)? shift @_: 0;

  my ($mean90a, $mean90p, $mean180, $varianceOrN90a, $varianceOrN90p, $varianceOrN180);

  if ($type eq "dev")
    {
    $varianceOrN90a = $self->numAngles() ;
    $varianceOrN90p = $self->numAngles() ;
    $varianceOrN180 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle90a) ** 2) / $varianceOrN90a ;} ($self->calcAllAngles90a($combo)) ; #ideal mean, overal variance
    map { $angleTestStat += (($_ - $expectedAngle90p) ** 2) / $varianceOrN90p ;} ($self->calcAllAngles90p($combo)) ;
    map { $angleTestStat += (($_ - $expectedAngle180) ** 2) / $varianceOrN180 ;} ($self->calcAllAngles180($combo)) ;

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean90a = ($$angleStats{"octahedral"}{"90"}{"mean"})? ($$angleStats{"octahedral"}{"90"}{"mean"}) : $expectedAngle90a ;
      $mean90p = ($$angleStats{"octahedral"}{"90"}{"mean"})? ($$angleStats{"octahedral"}{"90"}{"mean"}) : $expectedAngle90p ;
      #$varianceOrN90a = ($$angleStats{"octahedral"}{"90"}{"variance"})? ($$angleStats{"octahedral"}{"90"}{"variance"}) : $$angleStats{"variance"}; ## major variance
      #$varianceOrN90p = ($$angleStats{"octahedral"}{"90"}{"variance"})? ($$angleStats{"octahedral"}{"90"}{"variance"}) : $$angleStats{"variance"};
      $varianceOrN90a = $$angleStats{"variance"}; ## overall variance
      $varianceOrN90p = $$angleStats{"variance"};
      }
    elsif ($type eq "ownStats")
      {
      $mean90a = $$angleStats{"squarePyramidalVacancy"}{"90a"}{"mean"};
      $mean90p = $$angleStats{"squarePyramidalVacancy"}{"90p"}{"mean"};
      $varianceOrN90a = $$angleStats{"squarePyramidalVacancy"}{"90a"}{"variance"};
      $varianceOrN90p = $$angleStats{"squarePyramidalVacancy"}{"90p"}{"variance"};
      }

    #my @means = ($expectedAngle90a, $expectedAngle90p, $expectedAngle90p, $expectedAngle90p, $expectedAngle90p); ## ideal mean
    my @means = ($mean90a, $mean90p, $mean90p, $mean90p, $mean90p); ## major mean
    my @angles = ($self->calcAllAngles90a($combo), $self->calcAllAngles90p($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    ## set the smallest angle's diff to zero
    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std90a = 1/sqrt($varianceOrN90a);
    my $std90p = 1/sqrt($varianceOrN90p);
    my $invStds = [$std90a, $std90p, $std90p, $std90p, $std90p];
    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

    return $chiStat;
    }
  }

sub angleList
  {
  my $self = shift @_;

  return 0 if (! exists $self->{bestCombo});

  my @anglelist;
  push (@anglelist, $self->calcAllAngles90a());
  push (@anglelist, $self->calcAllAngles90p());
  push (@anglelist, $self->calcAllAngles180());

  my $pair1 = [sort {$a <=> $b} ($anglelist[0], $anglelist[5])];
  my $pair2 = [sort {$a <=> $b} ($anglelist[1], $anglelist[4])];
  my $pair3 = [sort {$a <=> $b} ($anglelist[2], $anglelist[3])];

  my @pairlist = sort { $$a[0] <=> $$b[0]; } ($pair1, $pair2, $pair3);
  @anglelist = ( @{$pairlist[0]}, @{$pairlist[1]}, @{$pairlist[2]} );

  return @anglelist;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"squarePyramidalVacancy"}{"90a"}}, $_->calcAllAngles90a());
  push (@{$$angleStats{"squarePyramidalVacancy"}{"90p"}}, $_->calcAllAngles90p());
  push (@{$$angleStats{"squarePyramidalVacancy"}{"180"}}, $_->calcAllAngles180());

  return 0;
  }


sub calcAllAngles90p
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};
  my $pair1 = [$$combo[0], $$combo[1]];
  my $pair2 = [$$combo[2], $$combo[3]];

  # angles between pairs
  for(my $x = 0; $x < @$pair1; $x++) 
    {
    for(my $y = 0; $y < @$pair2; $y++)
      {
      push (@angles, $center->angle($$pair1[$x], $$pair2[$y])) ;
      } 
    }
  return @angles;
  }


sub calcAllAngles90a
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};
      
  my @angles;
  my $center = $self->{shellObj}->{center};
  my $pair1 = [$$combo[0], $$combo[1]];

  push (@angles, $center->angle($$combo[0], $$combo[1])) ;

  return @angles;
  }


sub calcAllAngles180
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};
      
  my @angles;
  my $center = $self->{shellObj}->{center};
  my $pair2 = [$$combo[2], $$combo[3]];

  push (@angles, $center->angle($$combo[2], $$combo[3])) ;

  return @angles;
  }


sub calcAllAngles
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

