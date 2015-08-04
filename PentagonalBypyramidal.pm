## PentagonalBipyramidal.pm 

package PentagonalBypyramidal;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 7
			  );

our $expectedAngle72 = 72;
our $expectedAngle90 = 90;
our $expectedAngle144 = 144;
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
