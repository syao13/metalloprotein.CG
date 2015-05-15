## SquarePyramidal.pm

package SquarePyramidal;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
			  "numAtoms" => 5
			  );

our $expectedAngle90a = 90;
our $expectedAngle90p = 90;
our $expectedAngle180 = 180;

our $invCorrM = [
	[1.3333333,	0.6666667,	0.0000000,	0.0000000,	0.000,	0.000,	0.000,	0.000],
	[0.6666667,	1.3333333,	0.0000000,	0.0000000,	0.000,	0.000,	0.000,	0.000],
	[0.0000000,	0.0000000,	1.3333333,	0.6666667,	0.000,	0.000,	0.000,	0.000],
	[0.0000000,	0.0000000,	0.6666667,	1.3333333,	0.000,	0.000,	0.000,	0.000],
	[0.0000000,	0.0000000,	0.0000000,	0.0000000,	0.625,	-0.125,	-0.125,	-0.375],
	[0.0000000,	0.0000000,	0.0000000,	0.0000000,	-0.125,	0.625,	-0.375,	-0.125],
	[0.0000000,	0.0000000,	0.0000000,	0.0000000,	-0.125,	-0.375,	0.625,	-0.125],
	[0.0000000,	0.0000000,	0.0000000,	0.0000000,	-0.375,	-0.125,	-0.125,	0.625]
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

  foreach my $axialIndex (0..$#$combo)
    {
    my @planarIndeces = grep { $axialIndex != $_; } (0..$#$combo);
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,1,2,3)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,2,1,3)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,3,1,2)] );
    }

  return $orderedCombos;
  }

sub angleTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $combo = shift @_;
  my $angleStats = shift @_;

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
      $mean90a = ($$angleStats{"octahedral"}{"90"}{"mean"})? ($$angleStats{"octahedral"}{"90"}{"mean"}) : $expectedAngle90a;
      $mean90p = ($$angleStats{"octahedral"}{"90"}{"mean"})? ($$angleStats{"octahedral"}{"90"}{"mean"}) : $expectedAngle90p;
      #$varianceOrN90a = ($$angleStats{"octahedral"}{"90"}{"variance"})? ($$angleStats{"octahedral"}{"90"}{"variance"}) : $$angleStats{"variance"}; ## major variance
      #$varianceOrN90p = ($$angleStats{"octahedral"}{"90"}{"variance"})? ($$angleStats{"octahedral"}{"90"}{"variance"}) : $$angleStats{"variance"};
      $varianceOrN90a = $$angleStats{"variance"}; ## overall variance
      $varianceOrN90p = $$angleStats{"variance"};
      }
    else
      {
      $varianceOrN90a = $$angleStats{"squarePyramidal"}{"90a"}{"variance"};
      $varianceOrN90p = $$angleStats{"squarePyramidal"}{"90p"}{"variance"};
      }

    #my @means = ($expectedAngle90a, $expectedAngle90a, $expectedAngle90a, $expectedAngle90a, $expectedAngle90p, $expectedAngle90p, $expectedAngle90p, $expectedAngle90p); ## ideal mean
    my @means = ($mean90a, $mean90a, $mean90a, $mean90a, $mean90p, $mean90p, $mean90p, $mean90p); ## major mean
    my @angles = ($self->calcAllAngles90a($combo), $self->calcAllAngles90p($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    my $std90a = 1/sqrt($varianceOrN90a);
    my $std90p = 1/sqrt($varianceOrN90p);
    my $invStds = [$std90a, $std90a, $std90a, $std90a, $std90p, $std90p, $std90p, $std90p];

    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

    #print "mean 90, $expect90a, , $expect90p\nangles: ";
    #print map {"$_, "; } (@angles);
    #print "\nstd90, $std90a, $std90p\n";
    #print "covariance: ", $chiStat->element(1,1);#, "; previous: $angleTestStat\n";

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

  return @anglelist;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"squarePyramidal"}{"90a"}}, $_->calcAllAngles90a());
  push (@{$$angleStats{"squarePyramidal"}{"90p"}}, $_->calcAllAngles90p());
  push (@{$$angleStats{"squarePyramidal"}{"180"}}, $_->calcAllAngles180());

  return 0;
  }



sub calcAllAngles90p
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};
  my $axial = $$combo[0];
  my $pair1 = [$$combo[1], $$combo[2]];
  my $pair2 = [$$combo[3], $$combo[4]];

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
  my $axial = $$combo[0];
  my $pair1 = [$$combo[1], $$combo[2]];
  my $pair2 = [$$combo[3], $$combo[4]];

  # angles within a pair
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


sub calcAllAngles90a
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};
      
  my @angles;
  my $center = $self->{shellObj}->{center};
  my $axial = $$combo[0];
  my $pair1 = [$$combo[1], $$combo[2]];
  my $pair2 = [$$combo[3], $$combo[4]];

  # angles betwen axial and planar
  foreach my $pair ($pair1, $pair2)
    {
    for(my $x = 0; $x < @$pair; $x++)
      {
      push (@angles, $center->angle($axial, $$pair[$x])) ;
      }
    }

  return @angles;
  }


