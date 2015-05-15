## "RawStatistics.pm"
#
package RawStatistics;
use strict;


our @defaultDataMembers = (
                          "variables" => 0 # ref to a array of variable
                          );


sub new
  {
  my $class = shift @_;
  my $self = { @defaultDataMembers, @_ };
  bless $self, ref $class || $class;

  return $self;
  }


## Calculate mean
sub calcMean
  {
  my $self = shift @_;

  return 0 if ($self->{variables} == 0) ;
  my $variables = $self->{variables} ;

  my $mean;
  foreach my $value (@$variables)
    {
    $mean += $value;
    }
  
  my $n =  $self->count();
  $mean = $mean / $n;
  return $mean;
  }

## Calculate deviation
sub calcDeviation
  {
  my $self = shift @_;

  return 0 if ($self->{variables} == 0) ;
  my $variables = $self->{variables} ;

  my $mean = $self->calcMean();

  my $dev = 0;
  foreach my $value (@$variables)
    { $dev += ($mean - $value) ** 2; }

  return $dev;
  }


## Calculate variance
sub calcVariance
  {
  my $self = shift @_;

  return 0 if ($self->{variables} == 0) ;

  my $dev = $self->calcDeviation();
  my $n = $self->count();
  my $var = $dev / $n;
 
  return $var;
  }


## calculate standard deviation
sub calcStd
  {
  my $self = shift @_;

  return 0 if ($self->{variables} == 0) ;

  my $dev = $self->calcDeviation();
  my $n = $self->count();
  my $std = ($dev / $n) ** 0.5;

  return $std;
  }


## Calculate number of variables
sub count
  {
  my $self = shift @_;

  return 0 if ($self->{variables} == 0) ;
  my $variables = $self->{variables} ;

  return scalar(@$variables);
  }


## Find the minimum value
sub min
  {
  my $self = shift @_;

  return 0 if ($self->{variables} == 0) ;
  my @variables = @{$self->{variables}} ;

  my $min = (sort {$a <=> $b} (@variables))[0];

  return $min;
  }

## Find the maximum value
sub max
  {
  my $self = shift @_;

  return 0 if ($self->{variables} == 0) ;
  my @variables = @{$self->{variables}} ;

  my $max = (sort {$b <=> $a} (@variables))[0];

  return $max;
  }












 
