## SquareAntiprismaticMonocapped.pmm 

package SquareAntiprismaticMonocapped;
use strict;
use AtomShell;

use base "Coordination";

our @defaultDataMembers = (
                          "numAtoms" => 9
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
