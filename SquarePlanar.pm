## SquarePlanar.pm

package SquarePlanar;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 4
			  );

our $expectedAngle90 = 90;
our $expectedAngle180 = 180;

our $invCorrM = [
	[ 0.625, -0.125, -0.125, -0.375,    0,    0],
	[-0.125,  0.625, -0.375, -0.125,    0,    0],
	[-0.125, -0.375,  0.625, -0.125,    0,    0],
	[-0.375, -0.125, -0.125,  0.625,    0,    0],
	[0.000,  0.000,  0.000,  0.000,    1,    0],
	[0.000,  0.000,  0.000,  0.000,    0,    1]
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
   
  my $orderedCombo;
  push (@$orderedCombo, [$$combo[0], $$combo[1], $$combo[2], $$combo[3]]);
  push (@$orderedCombo, [$$combo[0], $$combo[2], $$combo[1], $$combo[3]]);
  push (@$orderedCombo, [$$combo[0], $$combo[3], $$combo[1], $$combo[2]]);

  return $orderedCombo;
  }

sub angleTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $combo = shift @_;
  my $angleStats = shift @_;
  my $leaveOut = (@_)? shift @_: 0;

  my ($mean90, $mean180, $varianceOrN90, $varianceOrN180);

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
      $mean180 = ($$angleStats{"octahedral"}{"180"}{"mean"})? ($$angleStats{"octahedral"}{"180"}{"mean"}) : $expectedAngle180;
      #$varianceOrN90 = ($$angleStats{"octahedral"}{"90"}{"variance"})? ($$angleStats{"octahedral"}{"90"}{"variance"}) : $$angleStats{"variance"}; ## major coordination variance
      #$varianceOrN180 = ($$angleStats{"octahedral"}{"90"}{"variance"})? ($$angleStats{"octahedral"}{"90"}{"variance"}) : $$angleStats{"variance"};
      $varianceOrN90 = $$angleStats{"variance"}; ## ovarall variance
      $varianceOrN180 = $$angleStats{"variance"};
      }
    elsif ($type eq "ownStats")
      {
      $mean90 =  $$angleStats{"squarePlanar"}{"90"}{"mean"};
      $mean180 = $$angleStats{"squarePlanar"}{"180"}{"mean"};
      $varianceOrN90 =  $$angleStats{"squarePlanar"}{"90"}{"variance"};
      $varianceOrN180 = $$angleStats{"squarePlanar"}{"180"}{"variance"};
      }

    #my @means = ($expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle90, $expectedAngle180, $expectedAngle180); ## ideal mean
    my @means = ($mean90, $mean90, $mean90, $mean90, $mean180, $mean180); ## major coordination mean
    my @angles = ($self->calcAllAngles90($combo), $self->calcAllAngles180($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    ## set the smallest angle's diff to zero
    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std90 = 1/sqrt($varianceOrN90);
    my $std180 = 1/sqrt($varianceOrN180);
    my $invStds = [$std90, $std90, $std90, $std90, $std180, $std180];
    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

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

  my $pair1 = [sort {$a <=> $b} ($anglelist[0], $anglelist[3])];
  my $pair2 = [sort {$a <=> $b} ($anglelist[1], $anglelist[2])];
  my $pair3 = [sort {$a <=> $b} ($anglelist[4], $anglelist[5])];

  my @pairlist = sort { $$a[0] <=> $$b[0]; } ($pair1, $pair2, $pair3);
  @anglelist = ( @{$pairlist[0]}, @{$pairlist[1]}, @{$pairlist[2]} );

  return @anglelist;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"squarePlanar"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"squarePlanar"}{"180"}}, $_->calcAllAngles180());

  return 0;
  }


sub calcAllAngles90
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


sub calcAllAngles180
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};
  my $pair1 = [$$combo[0], $$combo[1]];
  my $pair2 = [$$combo[2], $$combo[3]];

  # angles within pairs
  foreach my $pair ($pair1, $pair2)
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


