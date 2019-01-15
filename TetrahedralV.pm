## TetrahedralV.pm

############################################################################################
###
###   Written by Sen Yao, 07/20/2016
###   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/20/2016. All rights reserved.
###
############################################################################################

package TetrahedralV;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 3
			  );

our $expectedAngle = 109.5;

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
      $mean109 = ($$angleStats{"tetrahedral"}{"109"}{"mean"})? ($$angleStats{"tetrahedral"}{"109"}{"mean"}) : $expectedAngle ;
      #$varianceOrN = ($$angleStats{"tetrahedral"}{"109"}{"variance"})? ($$angleStats{"tetrahedral"}{"109"}{"variance"}) : $$angleStats{"variance"}; ## major variance
      $varianceOrN = $$angleStats{"variance"}; ## overall variance
      }
    else
      {
      $varianceOrN = $$angleStats{"tetrahedralVacancy"}{"109"}{"variance"};
      }

    my $chiStat;
    #map { $chiStat += (($_ - $expectedAngle) ** 2) / $varianceOrN ;} ($self->calcAllAngles109($combo)) ; ## ideal mean
    map { $chiStat += (($_ - $mean109) ** 2) / $varianceOrN ;} ($self->calcAllAngles109($combo)) ; ## major mean

    return $chiStat;
    }
  }

sub angleList
  {
  my $self = shift @_;

  return 0 if (! exists $self->{bestCombo});

  my @anglelist;
  push (@anglelist, $self->calcAllAngles109());

  return @anglelist;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"tetrahedralVacancy"}{"109"}}, $_->calcAllAngles109());

  return 0;
  }

sub calcAllAngles109
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  for(my $x = 0; $x < $#$combo; $x++)
    {
    for(my $y = $x+1; $y < @$combo; $y++)
      {
      push (@angles, $center->angle($$combo[$x], $$combo[$y])) ;
      } 
    }

  return @angles;
  }


