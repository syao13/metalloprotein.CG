## Cube.pm 

package Cube;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 8
			  );

our $expectedAngle70 = 70.5;
our $expectedAngle109 = 109.5;
our $expectedAngle180 = 180;


our $invCorrM = [
		];

sub new
  {
  my $class = shift @_;
  my $self = Coordination->new(@defaultDataMembers, @_);

  bless $self, ref $class || $class;

  return $self;
  }
