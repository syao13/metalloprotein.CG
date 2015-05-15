## TrigonalPlanar.pm

package TrigonalPlanar;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 3
			  );

our $expectedAngle = 120;

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

  my ($mean120, $varianceOrN);

  if ($type eq "dev")
    {
    $varianceOrN = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle) ** 2) / $varianceOrN ;} ($self->calcAllAngles120($combo)) ;

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean120 = ($$angleStats{"trigonalPlanar"}{"120"}{"mean"})? ($$angleStats{"trigonalPlanar"}{"120"}{"mean"}) : $expectedAngle  ;
      #$varianceOrN = ($$angleStats{"trigonalPlanar"}{"120"}{"variance"})? ($$angleStats{"trigonalPlanar"}{"120"}{"variance"}) : $$angleStats{"variance"}; ## major variance
      $varianceOrN = $$angleStats{"variance"}; ##overall variance
      }
    else
      {
      $varianceOrN = $$angleStats{"trigonalPlanar"}{"120"}{"variance"};
      }

    my $chiStat;
    #map { $chiStat += (($_ - $expectedAngle) ** 2) / $varianceOrN ;} ($self->calcAllAngles120($combo)) ; ## ideal mean
    map { $chiStat += (($_ - $mean120) ** 2) / $varianceOrN ;} ($self->calcAllAngles120($combo)) ; ## major mean

    return $chiStat;
    }
  }

sub angleList
  {
  my $self = shift @_;

  return 0 if (! exists $self->{bestCombo});

  my @anglelist;
  push (@anglelist, $self->calcAllAngles120());

  return @anglelist;
  }



sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"trigonalPlanar"}{"120"}}, $_->calcAllAngles120());

  return 0;
  }

sub calcAllAngles120
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


