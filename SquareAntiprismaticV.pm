## SquareAntiprismaticV.pm 

package SquareAntiprismaticV;
use strict;
use AtomShell;
use Time::HiRes qw(time);
use POSIX qw(strftime);

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 7,
                          "degreeFreedom" => 28
			  );

our $expectedAngle70 = 70.5;
our $expectedAngle82 = 82;
our $expectedAngle109 = 109.5;
our $expectedAngle143 = 143.6;


our $invCorrM = [
	[ 0.304, -0.054, -0.047, -0.052, -0.007, -0.008, -0.113, -0.037,  0.067, -0.039, -0.110,  0.064,  0.076,  0.079, -0.091,  0.102,  0.004,  0.011,  0.102,  0.011,  0.002],
	[-0.054,  0.311,  0.006, -0.063, -0.055, -0.047,  0.060, -0.045, -0.107,  0.110, -0.019,  0.048,  0.161,  0.040, -0.029,  0.009, -0.148,  0.013, -0.092,  0.076,  0.030],
	[-0.047,  0.006,  0.219,  0.004, -0.041, -0.041,  0.060,  0.053, -0.102,  0.056,  0.059, -0.100,  0.090,  0.088, -0.073,  0.006,  0.036, -0.073,  0.004, -0.072,  0.039],
	[-0.052, -0.063,  0.004,  0.312, -0.048, -0.055, -0.020,  0.108,  0.048, -0.045,  0.060, -0.112,  0.037,  0.157, -0.031, -0.091,  0.034,  0.077,  0.009,  0.013, -0.146],
	[-0.007, -0.055, -0.041, -0.048,  0.263,  0.038, -0.022, -0.165,  0.051, -0.118,  0.060, -0.058, -0.048, -0.102,  0.049,  0.028, -0.076,  0.023,  0.073, -0.108,  0.018],
	[-0.008, -0.047, -0.041, -0.055,  0.038,  0.261,  0.059, -0.117, -0.059, -0.165, -0.020,  0.050, -0.102, -0.045,  0.048,  0.072,  0.017, -0.108,  0.028,  0.024, -0.076],
	[-0.113,  0.060,  0.060, -0.020, -0.022,  0.059,  0.377, -0.036, -0.136,  0.058, -0.094,  0.095, -0.060,  0.051, -0.123, -0.051, -0.063, -0.030,  0.129, -0.012, -0.012],
	[-0.037, -0.045,  0.053,  0.108, -0.165, -0.117, -0.036,  0.390, -0.059,  0.191,  0.056, -0.107,  0.092, -0.021,  0.056, -0.044, -0.040,  0.033,  0.055,  0.108, -0.040],
	[ 0.067, -0.107, -0.102,  0.048,  0.051, -0.059, -0.136, -0.059,  0.330, -0.109,  0.095, -0.084, -0.152,  0.055, -0.083,  0.013,  0.044, -0.028, -0.025,  0.048,  0.047],
	[-0.039,  0.110,  0.056, -0.045, -0.118, -0.165,  0.058,  0.191, -0.109,  0.395, -0.037, -0.059, -0.017,  0.092,  0.058,  0.057, -0.043,  0.108, -0.046,  0.034, -0.040],
	[-0.110, -0.019,  0.059,  0.060,  0.060, -0.020, -0.094,  0.056,  0.095, -0.037,  0.375, -0.135,  0.047, -0.065, -0.123,  0.131, -0.010, -0.014, -0.049, -0.031, -0.060],
	[ 0.064,  0.048, -0.100, -0.112, -0.058,  0.050,  0.095, -0.107, -0.084, -0.059, -0.135,  0.331,  0.059, -0.147, -0.081, -0.028,  0.044,  0.050,  0.012, -0.028,  0.040],
	[ 0.076,  0.161,  0.090,  0.037, -0.048, -0.102, -0.060,  0.092, -0.152, -0.017,  0.047,  0.059,  0.416, -0.050, -0.112, -0.022, -0.114,  0.041,  0.045, -0.012,  0.038],
	[ 0.079,  0.040,  0.088,  0.157, -0.102, -0.045,  0.051, -0.021,  0.055,  0.092, -0.065, -0.147, -0.050,  0.407, -0.113,  0.047,  0.040, -0.013, -0.021,  0.039, -0.110],
	[-0.091, -0.029, -0.073, -0.031,  0.049,  0.048, -0.123,  0.056, -0.083,  0.058, -0.123, -0.081, -0.112, -0.113,  0.383, -0.089, -0.004,  0.028, -0.090,  0.029, -0.001],
	[ 0.102,  0.009,  0.006, -0.091,  0.028,  0.072, -0.051, -0.044,  0.013,  0.057,  0.131, -0.028, -0.022,  0.047, -0.089,  0.256, -0.011, -0.008,  0.003, -0.007, -0.106],
	[ 0.004, -0.148,  0.036,  0.034, -0.076,  0.017, -0.063, -0.040,  0.044, -0.043, -0.010,  0.044, -0.114,  0.040, -0.004, -0.011,  0.235, -0.007, -0.105, -0.082, -0.013],
	[ 0.011,  0.013, -0.073,  0.077,  0.023, -0.108, -0.030,  0.033, -0.028,  0.108, -0.014,  0.050,  0.041, -0.013,  0.028, -0.008, -0.007,  0.185, -0.007, -0.055, -0.083],
	[ 0.102, -0.092,  0.004,  0.009,  0.073,  0.028,  0.129,  0.055, -0.025, -0.046, -0.049,  0.012,  0.045, -0.021, -0.090,  0.003, -0.105, -0.007,  0.255, -0.007, -0.011],
	[ 0.011,  0.076, -0.072,  0.013, -0.108,  0.024, -0.012,  0.108,  0.048,  0.034, -0.031, -0.028, -0.012,  0.039,  0.029, -0.007, -0.082, -0.055, -0.007,  0.183, -0.007],
	[ 0.002,  0.030,  0.039, -0.146,  0.018, -0.076, -0.012, -0.040,  0.047, -0.040, -0.060,  0.040,  0.038, -0.110, -0.001, -0.106, -0.013, -0.083, -0.011, -0.007,  0.233]
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
  foreach my $topAtoms (&Coordination::_restrictedCombinations(4, @$combo))
    {
    my $orderedTop;
    foreach my $i (0..3) 
      {
      foreach my $j (0..3) 
	{
	next if ($j == $i);
	foreach my $k (0..3)
	  {
	  next if ($k == $i || $k ==$j);
	  my $idx = (grep {$_ != $i && $_ != $j && $_ != $k; } (0..3))[0];
	  push (@$orderedTop, [$$topAtoms[$i], $$topAtoms[$j], $$topAtoms[$k], $$topAtoms[$idx]]);
	  }
	}
      }
    my @bottomAtoms = grep { my $test = $_; ! (grep { $_ == $test; } (@$topAtoms)); } (@$combo);

    foreach my $top (@$orderedTop)
      {
      push (@$orderedCombos, [@$top, $bottomAtoms[0], $bottomAtoms[1], $bottomAtoms[2]]);
      push (@$orderedCombos, [@$top, $bottomAtoms[1], $bottomAtoms[2], $bottomAtoms[0]]);
      push (@$orderedCombos, [@$top, $bottomAtoms[2], $bottomAtoms[0], $bottomAtoms[1]]);
      }
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

  my ($mean70, $mean109, $mean82, $mean143, $varianceOrN70, $varianceOrN82, $varianceOrN109, $varianceOrN143);

  if ($type eq "dev")
    {
    $varianceOrN70 = $self->numAngles() ;
    $varianceOrN82 = $self->numAngles() ;
    $varianceOrN109 = $self->numAngles() ;
    $varianceOrN143 = $self->numAngles() ;

    my $angleTestStat;
    map { $angleTestStat += (($_ - $expectedAngle70) ** 2) / $varianceOrN70 ;} ($self->calcAllAngles70($combo)) ; #ideal mean, number of variable
    map { $angleTestStat += (($_ - $expectedAngle82) ** 2) / $varianceOrN82 ;} ($self->calcAllAngles82($combo)) ; 
    map { $angleTestStat += (($_ - $expectedAngle109) ** 2) / $varianceOrN109 ;} ($self->calcAllAngles109($combo)) ; 
    map { $angleTestStat += (($_ - $expectedAngle143) ** 2) / $varianceOrN143 ;} ($self->calcAllAngles143($combo)) ;

    return $angleTestStat;
    }
  else
    {
    if ($type eq "chi")
      {
      $mean70 = ($$angleStats{"SquareAntiprismatic"}{"70"}{"mean"})? ($$angleStats{"SquareAntiprismatic"}{"70"}{"mean"}) : $expectedAngle70; ## major coordinaiton mean
      $mean82 = ($$angleStats{"SquareAntiprismatic"}{"82"}{"mean"})? ($$angleStats{"SquareAntiprismatic"}{"82"}{"mean"}) : $expectedAngle82;
      $mean109 = ($$angleStats{"SquareAntiprismatic"}{"109"}{"mean"})? ($$angleStats{"SquareAntiprismatic"}{"109"}{"mean"}) : $expectedAngle109;
      $mean143 = ($$angleStats{"SquareAntiprismatic"}{"143"}{"mean"})? ($$angleStats{"SquareAntiprismatic"}{"143"}{"mean"}) : $expectedAngle143;
      $varianceOrN70 = $$angleStats{"variance"}; ## overall variance
      $varianceOrN82 = $$angleStats{"variance"}; 
      $varianceOrN109 = $$angleStats{"variance"}; 
      $varianceOrN143 = $$angleStats{"variance"};
      }
    else
      {
      $varianceOrN70 =  $$angleStats{"SquareAntiprismaticVacancy"}{"70"}{"variance"};
      $varianceOrN82 =  $$angleStats{"SquareAntiprismaticVacancy"}{"82"}{"variance"};
      $varianceOrN109 =  $$angleStats{"SquareAntiprismaticVacancy"}{"109"}{"variance"};
      $varianceOrN143 =  $$angleStats{"SquareAntiprismaticVacancy"}{"143"}{"variance"};
      }

    my @means = ($mean70, $mean70, $mean70, $mean70, $mean70, $mean70, $mean82, $mean82, $mean82, $mean82, $mean82, $mean82, $mean109, $mean109, $mean109, $mean143, $mean143, $mean143, $mean143, $mean143, $mean143); 
    my @angles = ($self->calcAllAngles70($combo), $self->calcAllAngles82($combo), $self->calcAllAngles109($combo), $self->calcAllAngles143($combo),);
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    if ($leaveOut eq "l" || $leaveOut eq "leave" || $leaveOut eq "leaveOut")
      {
      my $n = 0;
      for(my $x = 1; $x <= $#angles; $x++)
        { if ($angles[$x] < $angles[$n]) {$n = $x;} }
      $$diff[$n] = 0;
      }

    my $std70 = 1/sqrt($varianceOrN70);
    my $std82 = 1/sqrt($varianceOrN82);
    my $std109 = 1/sqrt($varianceOrN109);
    my $std143 = 1/sqrt($varianceOrN143);
    my $invStds = [$std70, $std70, $std70, $std70, $std70, $std70, $std82, $std82, $std82, $std82, $std82, $std82, $std109, $std109, $std109, $std143, $std143, $std143, $std143, $std143, $std143];
    #my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

    my $chiStat = 0;
    map {$chiStat += $$diff[$_] ** 2 / $$invStds[$_];} (0..(@$diff-1));

    #print "mean 82, $expect82\nangles: ";
    ##print map {"$_, "; } (@angles);
    ##print "\nstd82, $std82\n";
    ##print "covariance: ", $chiStat;#, "; previous: $angleTestStat\n";
    return $chiStat;
    }
  }


sub calcAllAngles70
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};
  push (@angles, $center->angle($$combo[0], $$combo[1]));
  push (@angles, $center->angle($$combo[1], $$combo[2]));
  push (@angles, $center->angle($$combo[2], $$combo[3]));
  push (@angles, $center->angle($$combo[3], $$combo[0]));
  push (@angles, $center->angle($$combo[4], $$combo[5]));
  push (@angles, $center->angle($$combo[5], $$combo[6]));

  return @angles;
  }

sub calcAllAngles82
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[0], $$combo[4]));
  push (@angles, $center->angle($$combo[1], $$combo[5]));
  push (@angles, $center->angle($$combo[2], $$combo[6]));
  push (@angles, $center->angle($$combo[0], $$combo[5]));
  push (@angles, $center->angle($$combo[1], $$combo[6]));
  push (@angles, $center->angle($$combo[3], $$combo[4]));

  return @angles;
  }

sub calcAllAngles109
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[0], $$combo[2]));
  push (@angles, $center->angle($$combo[1], $$combo[3]));
  push (@angles, $center->angle($$combo[4], $$combo[6]));

  return @angles;
  }

sub calcAllAngles143
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[0], $$combo[6]));
  push (@angles, $center->angle($$combo[2], $$combo[4]));
  push (@angles, $center->angle($$combo[3], $$combo[5]));
  push (@angles, $center->angle($$combo[1], $$combo[4]));
  push (@angles, $center->angle($$combo[2], $$combo[5]));
  push (@angles, $center->angle($$combo[3], $$combo[6]));

  return @angles;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"SquareAntiprismaticVacancy"}{"70"}}, $_->calcAllAngles70());
  push (@{$$angleStats{"SquareAntiprismaticVacancy"}{"82"}}, $_->calcAllAngles82());
  push (@{$$angleStats{"SquareAntiprismaticVacancy"}{"109"}}, $_->calcAllAngles109());
  push (@{$$angleStats{"SquareAntiprismaticVacancy"}{"143"}}, $_->calcAllAngles143());

  return 0;
  }


