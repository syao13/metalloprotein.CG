## SquareAntiprismaticBicapped.pmm 

############################################################################################
###
###   Written by Sen Yao, 07/20/2016
###   Copyright Sen Yao, Robert Flight, and Hunter Moseley, 07/20/2016. All rights reserved.
###
############################################################################################

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
