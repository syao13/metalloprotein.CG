## SquareAntiprismaticBicapped.pmm 

package SquareAntiprismaticBicapped;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 10 
			  );


our $invCorrM = [
		];

sub new
  {
  my $class = shift @_;
  my $self = Coordination->new(@defaultDataMembers, @_);

  bless $self, ref $class || $class;

  return $self;
  }
