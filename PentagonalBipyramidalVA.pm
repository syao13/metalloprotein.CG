## PentagonalBipyramidal.pm 

package PentagonalBipyramidalVA;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 6
			  );

our $expectedAngle72 = 72;
our $expectedAngle90 = 90;
our $expectedAngle144 = 144;


our $invCorrM = [
	[ 0.16, -0.08,  0.00,  0.00, -0.08,  0.08, -0.08,  0.00, -0.08,  0.08,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.08,  0.16, -0.08,  0.00,  0.00,  0.08,  0.08, -0.08,  0.00, -0.08,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.00, -0.08,  0.16, -0.08,  0.00, -0.08,  0.08,  0.08, -0.08,  0.00,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.00,  0.00, -0.08,  0.16, -0.08,  0.00, -0.08,  0.08,  0.08, -0.08,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.08,  0.00,  0.00, -0.08,  0.16, -0.08,  0.00, -0.08,  0.08,  0.08,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.08,  0.08, -0.08,  0.00, -0.08,  0.16,  0.00, -0.08, -0.08,  0.00,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.08,  0.08,  0.08, -0.08,  0.00,  0.00,  0.16,  0.00, -0.08, -0.08,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.00, -0.08,  0.08,  0.08, -0.08, -0.08,  0.00,  0.16,  0.00, -0.08,  0.000,  0.000,  0.000,  0.000,  0.000],
	[-0.08,  0.00, -0.08,  0.08,  0.08, -0.08, -0.08,  0.00,  0.16,  0.00,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.08, -0.08,  0.00, -0.08,  0.08,  0.00, -0.08, -0.08,  0.00,  0.16,  0.000,  0.000,  0.000,  0.000,  0.000],
	[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  1.414, -0.163,  0.456,  0.456, -0.163],
	[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -0.163,  1.414, -0.163,  0.456,  0.456],
	[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.456, -0.163,  1.414, -0.163,  0.456],
	[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.456,  0.456, -0.163,  1.414, -0.163],
	[ 0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00, -0.163,  0.456,  0.456, -0.163,  1.414],
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
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,1,2,3,4)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,1,2,4,3)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,1,3,2,4)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,1,3,4,2)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,1,4,2,3)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,1,4,3,2)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,2,1,3,4)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,2,1,4,3)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,2,3,1,4)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,2,4,1,3)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,3,1,2,4)] );
    push (@$orderedCombos, [$$combo[$axialIndex], map { $$combo[$planarIndeces[$_]]; } (0,3,2,1,4)] );
    }

  return $orderedCombos;
  }

sub angleTestStatistic
  {
  my $self = shift @_;
  my $type = shift @_;
  my $combo = shift @_;
  my $angleStats = shift @_;
  my $leaveOut = (@_)? shift @_: 0;

  my ($mean72, $mean144, $mean90, $mean180, $varianceOrN72, $varianceOrN90, $varianceOrN144, $varianceOrN180);

  if ($type eq "dev")
    {
    $varianceOrN72 = $self->numAngles() ;
    $varianceOrN90 = $self->numAngles() ;
    $varianceOrN144 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle72) ** 2) / $varianceOrN72 ;} ($self->calcAllAngles72($combo)) ; #ideal mean, number of variable
    map { $angleTestStat += (($_ - $expectedAngle90) ** 2) / $varianceOrN90 ;} ($self->calcAllAngles90($combo)) ; 
    map { $angleTestStat += (($_ - $expectedAngle144) ** 2) / $varianceOrN144 ;} ($self->calcAllAngles144($combo)) ; 

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean72 = ($$angleStats{"pentagonalBipyramidal"}{"72"}{"mean"})? ($$angleStats{"pentagonalBipyramidal"}{"72"}{"mean"}) : $expectedAngle72; ## major coordinaiton mean
      $mean90 = ($$angleStats{"pentagonalBipyramidal"}{"90"}{"mean"})? ($$angleStats{"pentagonalBipyramidal"}{"90"}{"mean"}) : $expectedAngle90;
      $mean144 = ($$angleStats{"pentagonalBipyramidal"}{"144"}{"mean"})? ($$angleStats{"pentagonalBipyramidal"}{"144"}{"mean"}) : $expectedAngle144;
      $varianceOrN72 = $$angleStats{"variance"}; ## overall variance
      $varianceOrN90 = $$angleStats{"variance"}; 
      $varianceOrN144 = $$angleStats{"variance"}; 
      }
    else
      {
      $varianceOrN72 =  $$angleStats{"pentagonalBipyramidalVacancyAxial"}{"72"}{"variance"};
      $varianceOrN90 =  $$angleStats{"pentagonalBipyramidalVacancyAxial"}{"90"}{"variance"};
      $varianceOrN144 =  $$angleStats{"pentagonalBipyramidalVacancyAxial"}{"144"}{"variance"};
      }

    my @means = ($mean72, $mean72, $mean72, $mean72, $mean72, $mean90, $mean90, $mean90, $mean90, $mean90, $mean144, $mean144, $mean144, $mean144, $mean144); 
    my @angles = ($self->calcAllAngles72($combo), $self->calcAllAngles90($combo), $self->calcAllAngles144($combo));
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std72 = 1/sqrt($varianceOrN72);
    my $std90 = 1/sqrt($varianceOrN90);
    my $std144 = 1/sqrt($varianceOrN144);
    my $invStds = [$std72, $std72, $std72, $std72, $std72, $std90, $std90, $std90, $std90, $std90, $std144, $std144, $std144, $std144, $std144];

    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

    #print "mean 90, $expect90\nangles: ";
    ##print map {"$_, "; } (@angles);
    ##print "\nstd90, $std90\n";
    ##print "covariance: ", $chiStat;#, "; previous: $angleTestStat\n";
    
    return $chiStat;
    }
  }

sub calcAllAngles90
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0]];
  my $planar = [$$combo[1], $$combo[2], $$combo[3], $$combo[4], $$combo[5]];
  my $center = $self->{shellObj}->{center};

  for(my $x = 0; $x < @$axial; $x++) # axial atom to planar atom angles, 90 is ideal
    {
    for(my $y = 0; $y < @$planar; $y++)
      {
      push (@angles, $center->angle($$axial[$x], $$planar[$y]));
      }
    }

  return @angles;
  }

sub calcAllAngles72
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0]];
  my $planar = [$$combo[1], $$combo[2], $$combo[3], $$combo[4], $$combo[5]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$planar[0], $$planar[1])); ## planar atom to adjacent planar atom angles, 72 is ideal
  push (@angles, $center->angle($$planar[1], $$planar[2])); 
  push (@angles, $center->angle($$planar[2], $$planar[3])); 
  push (@angles, $center->angle($$planar[3], $$planar[4])); 
  push (@angles, $center->angle($$planar[4], $$planar[0])); 

  return @angles;
  }

sub calcAllAngles144
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $axial = [$$combo[0]];
  my $planar = [$$combo[1], $$combo[2], $$combo[3], $$combo[4], $$combo[5]];
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$planar[0], $$planar[2])); ## planar atom to every other planar atom angles, 144 is ideal
  push (@angles, $center->angle($$planar[1], $$planar[3])); 
  push (@angles, $center->angle($$planar[2], $$planar[4])); 
  push (@angles, $center->angle($$planar[3], $$planar[0])); 
  push (@angles, $center->angle($$planar[4], $$planar[1]));

  return @angles;
  }

sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"pentagonalBipyramidalVacancyAxial"}{"72"}}, $_->calcAllAngles72());
  push (@{$$angleStats{"pentagonalBipyramidalVacancyAxial"}{"90"}}, $_->calcAllAngles90());
  push (@{$$angleStats{"pentagonalBipyramidalVacancyAxial"}{"144"}}, $_->calcAllAngles144());

  return 0;
  }


