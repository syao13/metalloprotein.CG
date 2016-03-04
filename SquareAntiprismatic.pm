## SquareAntiprismatic.pm 

package SquareAntiprismatic;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 8
			  );

our $expectedAngle70 = 70.5;
our $expectedAngle82 = 82;
our $expectedAngle109 = 109.5;
our $expectedAngle143 = 143.6;


our $invCorrM = [
	[ 0.208, -0.017, -0.037, -0.017, -0.012, -0.013, -0.030, -0.029, -0.088, -0.040,  0.047,  0.031, -0.041, -0.087,  0.034,  0.046,  0.072,  0.071, -0.063, -0.019,  0.052, -0.081,  0.004,  0.011, -0.080,  0.053,  0.011,  0.003],
	[-0.017,  0.209, -0.017, -0.037, -0.029, -0.012, -0.011, -0.029,  0.033, -0.087, -0.040,  0.047,  0.048, -0.041, -0.088,  0.035,  0.070,  0.070, -0.018, -0.063,  0.011,  0.052, -0.082,  0.004,  0.004, -0.081,  0.054,  0.011],
	[-0.037, -0.017,  0.211, -0.017, -0.030, -0.029, -0.007, -0.011,  0.049,  0.034, -0.088, -0.040,  0.036,  0.049, -0.047, -0.087,  0.068,  0.068, -0.062, -0.019,  0.005,  0.011,  0.053, -0.081,  0.012,  0.003, -0.081,  0.053],
	[-0.017, -0.037, -0.017,  0.210, -0.012, -0.031, -0.029, -0.007, -0.040,  0.050,  0.033, -0.089, -0.086,  0.035,  0.047, -0.047,  0.070,  0.069, -0.021, -0.063, -0.081,  0.004,  0.012,  0.052,  0.052,  0.012,  0.002, -0.080],
	[-0.012, -0.029, -0.030, -0.012,  0.209, -0.017, -0.036, -0.016, -0.040, -0.087,  0.032,  0.048, -0.040,  0.048,  0.033, -0.088, -0.019, -0.063,  0.070,  0.071,  0.012,  0.003, -0.081,  0.051,  0.011,  0.053, -0.082,  0.005],
	[-0.013, -0.012, -0.029, -0.031, -0.017,  0.208, -0.017, -0.036,  0.048, -0.040, -0.088,  0.031, -0.088, -0.040,  0.048,  0.033, -0.061, -0.018,  0.070,  0.070,  0.052,  0.011,  0.003, -0.082,  0.004,  0.012,  0.052, -0.081],
	[-0.030, -0.011, -0.007, -0.029, -0.036, -0.017,  0.209, -0.017,  0.033,  0.049, -0.039, -0.090,  0.034, -0.087, -0.042,  0.049, -0.019, -0.063,  0.071,  0.070, -0.082,  0.052,  0.011,  0.005, -0.080,  0.003,  0.013,  0.051],
	[-0.029, -0.029, -0.011, -0.007, -0.016, -0.036, -0.017,  0.210, -0.089,  0.034,  0.049, -0.040,  0.049,  0.034, -0.088, -0.041, -0.064, -0.021,  0.071,  0.070,  0.005, -0.081,  0.052,  0.012,  0.052, -0.081,  0.004,  0.013],
	[-0.088,  0.033,  0.049, -0.040, -0.040,  0.048,  0.033, -0.089,  0.277, -0.025, -0.096, -0.027,  0.077, -0.067, -0.066,  0.080, -0.088,  0.042, -0.087,  0.038, -0.012, -0.020, -0.011, -0.022,  0.069,  0.067, -0.011, -0.016],
	[-0.040, -0.087,  0.034,  0.050, -0.087, -0.040,  0.049,  0.034, -0.025,  0.278, -0.024, -0.096,  0.080,  0.079, -0.067, -0.066,  0.038, -0.090,  0.040, -0.088, -0.023, -0.010, -0.021, -0.010, -0.014,  0.067,  0.068, -0.013],
	[ 0.047, -0.040, -0.088,  0.033,  0.032, -0.088, -0.039,  0.049, -0.096, -0.024,  0.276, -0.026, -0.064,  0.079,  0.077, -0.068, -0.086,  0.038, -0.089,  0.043, -0.010, -0.024, -0.010, -0.022, -0.015, -0.012,  0.066,  0.071],
	[ 0.031,  0.047, -0.040, -0.089,  0.048,  0.031, -0.090, -0.040, -0.027, -0.096, -0.026,  0.275, -0.067, -0.065,  0.079,  0.075,  0.042, -0.084,  0.041, -0.088, -0.022, -0.011, -0.024, -0.012,  0.069, -0.013, -0.014,  0.068],
	[-0.041,  0.048,  0.036, -0.086, -0.040, -0.088,  0.034,  0.049,  0.077,  0.080, -0.064, -0.067,  0.281, -0.024, -0.100, -0.026, -0.092,  0.039,  0.043, -0.088,  0.069, -0.014, -0.014,  0.068, -0.010, -0.026, -0.009, -0.020],
	[-0.087, -0.041,  0.049,  0.035,  0.048, -0.040, -0.087,  0.034, -0.067,  0.079,  0.079, -0.065, -0.024,  0.280, -0.026, -0.100,  0.039, -0.088, -0.089,  0.040,  0.070,  0.068, -0.012, -0.014, -0.021, -0.010, -0.026, -0.009],
	[ 0.034, -0.088, -0.047,  0.047,  0.033,  0.048, -0.042, -0.088, -0.066, -0.067,  0.077,  0.079, -0.100, -0.026,  0.281, -0.027, -0.085,  0.041,  0.039, -0.088, -0.014,  0.068,  0.069, -0.014, -0.011, -0.019, -0.013, -0.023],
	[ 0.046,  0.035, -0.087, -0.047, -0.088,  0.033,  0.049, -0.041,  0.080, -0.066, -0.068,  0.075, -0.026, -0.100, -0.027,  0.282,  0.041, -0.088, -0.086,  0.042, -0.016, -0.013,  0.066,  0.070, -0.023, -0.011, -0.018, -0.013],
	[ 0.072,  0.070,  0.068,  0.070, -0.019, -0.061, -0.019, -0.064, -0.088,  0.038, -0.086,  0.042, -0.092,  0.039, -0.085,  0.041,  0.315, -0.033, -0.083, -0.084, -0.040,  0.027, -0.039,  0.028, -0.038,  0.026, -0.039,  0.023],
	[ 0.071,  0.070,  0.068,  0.069, -0.063, -0.018, -0.063, -0.021,  0.042, -0.090,  0.038, -0.084,  0.039, -0.088,  0.041, -0.088, -0.033,  0.315, -0.084, -0.083,  0.027, -0.038,  0.027, -0.041,  0.025, -0.039,  0.024, -0.038],
	[-0.063, -0.018, -0.062, -0.021,  0.070,  0.070,  0.071,  0.071, -0.087,  0.040, -0.089,  0.041,  0.043, -0.089,  0.039, -0.086, -0.083, -0.084,  0.316, -0.030, -0.038,  0.025, -0.040,  0.025,  0.024, -0.039,  0.026, -0.037],
	[-0.019, -0.063, -0.019, -0.063,  0.071,  0.070,  0.070,  0.070,  0.038, -0.088,  0.043, -0.088, -0.088,  0.040, -0.088,  0.042, -0.084, -0.083, -0.030,  0.316,  0.026, -0.040,  0.024, -0.038, -0.039,  0.025, -0.037,  0.025],
	[ 0.052,  0.011,  0.005, -0.081,  0.012,  0.052, -0.082,  0.005, -0.012, -0.023, -0.010, -0.022,  0.069,  0.070, -0.014, -0.016, -0.040,  0.027, -0.038,  0.026,  0.169, -0.005, -0.007, -0.006, -0.068, -0.002, -0.002, -0.068],
	[-0.081,  0.052,  0.011,  0.004,  0.003,  0.011,  0.052, -0.081, -0.020, -0.010, -0.024, -0.011, -0.014,  0.068,  0.068, -0.013,  0.027, -0.038,  0.025, -0.040, -0.005,  0.169, -0.004, -0.007, -0.068, -0.068, -0.002, -0.002],
	[ 0.004, -0.082,  0.053,  0.012, -0.081,  0.003,  0.011,  0.052, -0.011, -0.021, -0.010, -0.024, -0.014, -0.012,  0.069,  0.066, -0.039,  0.027, -0.040,  0.024, -0.007, -0.004,  0.170, -0.005, -0.001, -0.069, -0.069, -0.002],
	[ 0.011,  0.004, -0.081,  0.052,  0.051, -0.082,  0.005,  0.012, -0.022, -0.010, -0.022, -0.012,  0.068, -0.014, -0.014,  0.070,  0.028, -0.041,  0.025, -0.038, -0.006, -0.007, -0.005,  0.170, -0.003, -0.001, -0.068, -0.069],
	[-0.080,  0.004,  0.012,  0.052,  0.011,  0.004, -0.080,  0.052,  0.069, -0.014, -0.015,  0.069, -0.010, -0.021, -0.011, -0.023, -0.038,  0.025,  0.024, -0.039, -0.068, -0.068, -0.001, -0.003,  0.168, -0.004, -0.008, -0.003],
	[ 0.053, -0.081,  0.003,  0.012,  0.053,  0.012,  0.003, -0.081,  0.067,  0.067, -0.012, -0.013, -0.026, -0.010, -0.019, -0.011,  0.026, -0.039, -0.039,  0.025, -0.002, -0.068, -0.069, -0.001, -0.004,  0.169, -0.004, -0.008],
	[ 0.011,  0.054, -0.081,  0.002, -0.082,  0.052,  0.013,  0.004, -0.011,  0.068,  0.066, -0.014, -0.009, -0.026, -0.013, -0.018, -0.039,  0.024,  0.026, -0.037, -0.002, -0.002, -0.069, -0.068, -0.008, -0.004,  0.169, -0.003],
	[ 0.003,  0.011,  0.053, -0.080,  0.005, -0.081,  0.051,  0.013, -0.016, -0.013,  0.071,  0.068, -0.020, -0.009, -0.023, -0.013,  0.023, -0.038, -0.037,  0.025, -0.068, -0.002, -0.002, -0.069, -0.003, -0.008, -0.003,  0.168]
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
  my %tops;
  foreach my $topAtoms (&Coordination::_restrictedCombinations(4, @$combo))
    {
    next if ($tops{join ("", sort (@$topAtoms))});

    my $orderedTop;
    push (@$orderedTop, [$$topAtoms[0], $$topAtoms[1], $$topAtoms[2], $$topAtoms[3]]);
    push (@$orderedTop, [$$topAtoms[0], $$topAtoms[2], $$topAtoms[3], $$topAtoms[1]]);
    push (@$orderedTop, [$$topAtoms[0], $$topAtoms[3], $$topAtoms[1], $$topAtoms[2]]);
    my @bottomAtoms = grep { my $test = $_; ! (grep { $_ == $test; } (@$topAtoms)); } (@$combo);

    my $bottomString = join ("", sort (@bottomAtoms));
    $tops{$bottomString} = 1;

    foreach my $top (@$orderedTop)
      {
      for my $i (0..3) 
	{
	for my $j (0..3)
	  {
	  next if ($i == $j);
	  for my $k (0..3) 
	    {
	    next if ($i == $k && $j == $k);
	      {
	      my $l = (grep {$_ != $i && $_ != $j && $_ != $k ;} (0..3))[0];
	      push (@$orderedCombos, [@$top, $bottomAtoms[$i], $bottomAtoms[$j], $bottomAtoms[$k], $bottomAtoms[$l]]);
	      }
	    }
	  }
	}
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
      $varianceOrN70 =  $$angleStats{"SquareAntiprismatic"}{"70"}{"variance"};
      $varianceOrN82 =  $$angleStats{"SquareAntiprismatic"}{"82"}{"variance"};
      $varianceOrN109 =  $$angleStats{"SquareAntiprismatic"}{"109"}{"variance"};
      $varianceOrN143 =  $$angleStats{"SquareAntiprismatic"}{"143"}{"variance"};
      }

    my @means = ($mean70, $mean70, $mean70, $mean70, $mean70, $mean70, $mean70, $mean70, $mean82, $mean82, $mean82, $mean82, $mean82, $mean82, $mean82, $mean82, $mean109, $mean109, $mean109, $mean109, $mean143, $mean143, $mean143, $mean143, $mean143, $mean143, $mean143, $mean143); 
    my @angles = ($self->calcAllAngles70($combo), $self->calcAllAngles82($combo), $self->calcAllAngles109($combo), $self->calcAllAngles143($combo),);
    my $diff = [map { $angles[$_] - $means[$_]; } (0..(@angles-1))];

    my $std70 = 1/sqrt($varianceOrN70);
    my $std82 = 1/sqrt($varianceOrN82);
    my $std109 = 1/sqrt($varianceOrN109);
    my $std143 = 1/sqrt($varianceOrN143);
    my $invStds = [$std70, $std70, $std70, $std70, $std70, $std70, $std70, $std70, $std82, $std82, $std82, $std82, $std82, $std82, $std82, $std82, $std109, $std109, $std109, $std109, $std143, $std143, $std143, $std143, $std143, $std143, $std143, $std143];

    my $chiStat = $self->covMatChi($diff, $invStds, $invCorrM);

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
  push (@angles, $center->angle($$combo[6], $$combo[7]));
  push (@angles, $center->angle($$combo[7], $$combo[4]));

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
  push (@angles, $center->angle($$combo[3], $$combo[7]));
  push (@angles, $center->angle($$combo[0], $$combo[5]));
  push (@angles, $center->angle($$combo[1], $$combo[6]));
  push (@angles, $center->angle($$combo[2], $$combo[7]));
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
  push (@angles, $center->angle($$combo[5], $$combo[7]));

  return @angles;
  }

sub calcAllAngles143
  {
  my $self = shift @_;
  my $combo = (@_)? shift @_ : $self->{bestCombo}->{ligands};

  my @angles;
  my $center = $self->{shellObj}->{center};

  push (@angles, $center->angle($$combo[0], $$combo[6]));
  push (@angles, $center->angle($$combo[1], $$combo[7]));
  push (@angles, $center->angle($$combo[2], $$combo[4]));
  push (@angles, $center->angle($$combo[3], $$combo[5]));
  push (@angles, $center->angle($$combo[0], $$combo[7]));
  push (@angles, $center->angle($$combo[1], $$combo[4]));
  push (@angles, $center->angle($$combo[2], $$combo[5]));
  push (@angles, $center->angle($$combo[3], $$combo[6]));

  return @angles;
  }


sub calcExpectedAngleStats
  {
  my $self = shift @_;
  my $angleStats = shift @_;

  push (@{$$angleStats{"SquareAntiprismatic"}{"70"}}, $_->calcAllAngles70());
  push (@{$$angleStats{"SquareAntiprismatic"}{"82"}}, $_->calcAllAngles82());
  push (@{$$angleStats{"SquareAntiprismatic"}{"109"}}, $_->calcAllAngles109());
  push (@{$$angleStats{"SquareAntiprismatic"}{"143"}}, $_->calcAllAngles143());

  return 0;
  }


