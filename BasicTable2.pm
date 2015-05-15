#       Copyright Hunter Moseley, 2003. All rights reserved.
#       Written by Hunter Moseley 2/23/2003
#	Modified by Hunter Moseley 11/16/2010
#
#  BasicTable.pm
#	Contains subroutines for reading and writing name-value pairs in a file format.
#
#	Subroutines:
#		readTableFile - read a table file and return a hash structure containing name-value pairs.
#		writeTableFile - write a table file using name-value pairs in a hash structure.
#		cloneTable - clones a table hash structure.
#		getTableValue - gets a specific value from the table hash structure.
#		getTableValueOrDefault - gets a specific value from the table hash if it exists or returns the default value.
#		testTableValue - tests if the key(s) exists in the table hash.
#		testTableStructure - tests the table hash against a template hash.
#
package BasicTable2;
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw();
@EXPORT_OK = qw(cloneTable readTableFile writeTableFile getTableValue getTableValueOrDefault testTableValue numTableValue testTableStructure);
%EXPORT_TAGS = ( ALL => [@EXPORT_OK] );

use strict;
#use Dumpvalue qw(:ALL);
#my $dumper    = new Dumpvalue;





# cloneTable
#   clones a table hash.
#
#        $original_hlist - the original $table_hlist
#        $cloned_hlist   - the cloned $table_hlist (optional).
sub cloneTable
  {
  my $original_hlist = shift @_;
  my $cloned_hlist = {};
  if (@_)
    { $cloned_hlist = shift @_; }

  # generic cloning code
  my @clone_stack;
  push @clone_stack, { "type" => "HASH", "clone" => $cloned_hlist, "range" => [ keys %$original_hlist ], "index" => 0, "original" => $original_hlist }; 

  while(@clone_stack)
    {
    my $index = $clone_stack[$#clone_stack]{"index"};

    if ($index < @{$clone_stack[$#clone_stack]{"range"}}) # clone a part of the current hash or array
      {
      # get new value
      my $test_value;
      if ($clone_stack[$#clone_stack]{"type"} eq "HASH")
	{ $test_value = $clone_stack[$#clone_stack]{"original"}{$clone_stack[$#clone_stack]{"range"}[$index]}; }
      elsif ($clone_stack[$#clone_stack]{"type"} eq "ARRAY")
	{ $test_value = $clone_stack[$#clone_stack]{"original"}[$index]; }

      my $test_ref = ref($test_value);
      if ($test_ref eq "HASH") # test value is a hash
	{ push @clone_stack, {"type" => "HASH", "clone" => {}, "range" => [ keys %$test_value ], "index" => 0, "original" => $test_value }; }
      elsif ($test_ref eq "ARRAY") # test value is an array
	{ push @clone_stack, {"type" => "ARRAY", "clone" => [], "range" => $test_value, "index" => 0, "original" => $test_value }; }
      elsif ($clone_stack[$#clone_stack]{"type"} eq "HASH") # clone section is a hash
	{
	$clone_stack[$#clone_stack]{"clone"}{$clone_stack[$#clone_stack]{"range"}[$index]} = $test_value;
	$clone_stack[$#clone_stack]{"index"}++;
	}
      else # clone section is an array
	{
	$clone_stack[$#clone_stack]{"clone"}[$index] = $test_value;
	$clone_stack[$#clone_stack]{"index"}++;
	}
      }
    else # pop finished cloned section.
      {
      if (@clone_stack > 1)
	{
	my $index = $clone_stack[$#clone_stack - 1]{"index"};

	if($clone_stack[$#clone_stack - 1]{"type"} eq "HASH") # handle hashes
	  { $clone_stack[$#clone_stack - 1]{"clone"}{$clone_stack[$#clone_stack - 1]{"range"}[$index]} = $clone_stack[$#clone_stack]{"clone"}; }
	else # handle arrays
	  { $clone_stack[$#clone_stack - 1]{"clone"}[$index] = $clone_stack[$#clone_stack]{"clone"}; }
	}
      
      pop @clone_stack;
      if (@clone_stack)
	{ $clone_stack[$#clone_stack]{"index"}++; }
      }
    }
  
  return $cloned_hlist;
  }


# readTableFile
#   reads a table file and returns a hash of name-value pairs.
#
#
#   Parameters:
#       $filename - table file to read.
#	$table_hash - ref to hash to put name-value pairs in (optional).
#
sub readTableFile
  {
  my $filename = shift @_;
  my $table_hash = {};
  if (@_)
    { $table_hash = shift @_; }

  local *TABLEFILE;
  if ($filename eq "-")
    { *TABLEFILE = *STDIN; }
  elsif ($filename =~ /[.]gz$/)
    { open (TABLEFILE, "zcat $filename|") || die "unable to open $filename"; }
  elsif ($filename =~ /[.]bz2$/)
    { open (TABLEFILE, "bzcat $filename|") || die "unable to open $filename"; }
  else
    { open (TABLEFILE, "<$filename") || die "unable to open $filename"; }

  my @hash_array;
  push @hash_array, [ "", $table_hash, 0 ];

  # divide TABLEFILE into save_frames.
  my $line_count = 0;
  while (my $line = <TABLEFILE>)
    {
    $line_count++;
    next if (($line =~ /^\s*\#/) || ($line =~ /^\s*$/)); # skip comments and blank lines.

    my @tokens = grep { $_ !~ /^\s*$/; } (split(/([\"][^\"]*[\"]|\s+|[:][:]?|[,\#\[\]])/,$line));

    my $name = shift @tokens;
    $name =~ s/^[\"](.*)[\"]$/$1/; # remove quotes

    if ($name eq "]") # handle end of object
      {
      my $object = pop @hash_array;
      if (! @hash_array)
	{ die "Improper \"]\" at line $line_count in $filename\n"; }
      
      &addPair($hash_array[$#hash_array][1],$$object[0], $$object[1], $$object[2]);
      next;
      }

    if ($tokens[0] !~ /^[:]/)
      { die "Missing \":\" at line $line_count in $filename\n"; }
    
    my $explicit_array = ($tokens[0] eq "::");
    shift @tokens; # skip colon

    while(@tokens)
      {
      my $single_token = shift @tokens;
      last if ($single_token eq "#"); # skip comments
      next if ($single_token =~ /^\s*$/); # skip tokens
      next if ($single_token eq ","); # skip commas

      if ($single_token eq "[") # create object
	{
	push @hash_array, [ $name , {}, $explicit_array ];
	last;
	}
      
      if ($single_token eq "]") # handle end of object
	{
	my $object = pop @hash_array;
	if (! @hash_array)
	  { die "Improper \"]\" at line $line_count in $filename\n"; }

	&addPair($hash_array[$#hash_array][1],$$object[0], $$object[1], $$object[2]);
	last;
	}
      
      $single_token =~ s/^[\"](.*)[\"]$/$1/; # remove quotes

      if ($single_token =~ /[\"]/)
	{ die "Improper double quote at line $line_count in $filename\n"; }

      # add the name value pair
      &addPair($hash_array[$#hash_array][1],$name,$single_token, $explicit_array);
      }
    
    }
  
  if (@hash_array > 1)
    { die "Missing \"]\" in $filename\n"; }

  close TABLEFILE if ($filename ne "-");

  return $table_hash; 
  }

# addPair
#   Add name value pair to hash
#
#   Parameters:
#	$hash - ref to hash.
#	$name - name of pair.
#	$value - value of pair.
#	$explicit_array - whether the value is part of an explicitly declared array.
#
sub addPair
  {
  my $hash = shift @_;
  my $name = shift @_;
  my $value = shift @_;
  my $explicit_array = shift @_;

  if (! exists $$hash{$name})
    {
    if (! $explicit_array)
      { $$hash{$name} = $value; }
    else
      { $$hash{$name} = [ $value ]; }

    return;
    }

  if (ref $$hash{$name} ne "ARRAY")
    { $$hash{$name} = [ $$hash{$name} ]; }

  push @{$$hash{$name}}, $value;
  return;
  }


# writeTableFile
#   Write out a table file from a given hash.
#
#   Parameters:
#	$filename - name of output filename.
#	$table_hlist - reference to bmrb hash structure.
#	$write_options - ref to hash of options (optional);
#
sub writeTableFile
  {
  my $filename = shift @_; 
  my $table_hash = shift @_;
  my $write_options = {};
  if (@_)
    { $write_options = shift @_; }
  

  my $order = [];
  my $order_hash = "";
  if (exists $$write_options{"print_order"} && (ref $$write_options{"print_order"} eq "ARRAY"))
    { 
    $order = $$write_options{"print_order"}[0]; 
    $order_hash = $$write_options{"print_order"}[1]; 
    }
  
  # start writing BMRB file
  local *TABLEFILE;
  my $open_success = 1;
  if ($filename eq "-")
    { *TABLEFILE = *STDOUT; }
  elsif ($filename =~ /[.]gz$/)
    { $open_success = open (TABLEFILE, "|gzip > $filename"); }
  elsif ($filename =~ /[.]bz2$/)
    { $open_success = open (TABLEFILE, "|bzip2 > $filename"); }
  else
    { $open_success = open (TABLEFILE, ">$filename"); }
  
  if (! $open_success && ! exists $$write_options{"no_die"})
    { die "unable to open $filename"; }
  elsif (! $open_success)
    { 
    print STDERR "unable to open $filename\n";
    return 0;
    }

  @$order = (grep { exists $$table_hash{$_}; } (@$order));

  if (! exists $$write_options{"limit_to_order"})
    { push @$order, (grep { my $test = $_; ! (grep {$_ eq $test; } (@$order)); } (sort {$a cmp $b; } (keys %$table_hash))); }

  my @print_stack;
  push @print_stack, { "type" => "HASH", "range" => [ @$order ], "index" => 0, "original" => $table_hash, "order_hash" => $order_hash }; 
  my $hash_depth = 0;

  while(@print_stack)
    {
    my $index = $print_stack[$#print_stack]{"index"};

    if ($index < @{$print_stack[$#print_stack]{"range"}}) # print a part of the current hash or array
      {
      # get new value
      my $test_value;
      if ($print_stack[$#print_stack]{"type"} eq "HASH")
	{ $test_value = $print_stack[$#print_stack]{"original"}{$print_stack[$#print_stack]{"range"}[$index]}; }
      elsif ($print_stack[$#print_stack]{"type"} eq "ARRAY")
	{ $test_value = $print_stack[$#print_stack]{"original"}[$index]; }

      my $test_ref = ref($test_value);
      if ($test_ref eq "HASH") # test value is a hash
	{ 
	if ($print_stack[$#print_stack]{"type"} eq "HASH")
	  { print TABLEFILE "  "x($hash_depth),$print_stack[$#print_stack]{"range"}[$index], " :"; }
	print TABLEFILE " [\n";
	if (exists $print_stack[$#print_stack]{"order_hash"} && (ref($print_stack[$#print_stack]{"order_hash"}) eq "HASH") && (exists $print_stack[$#print_stack]{"order_hash"}{$print_stack[$#print_stack]{"range"}[$index]} || exists $print_stack[$#print_stack]{"order_hash"}{$print_stack[$#print_stack]{"array_name"}}))
	  {	  
	  my $order, $order_hash;
	  if ($print_stack[$#print_stack]{"type"} eq "HASH")
	    {
	    $order = $print_stack[$#print_stack]{"order_hash"}{$print_stack[$#print_stack]{"range"}[$index]}[0];
	    $order_hash = $print_stack[$#print_stack]{"order_hash"}{$print_stack[$#print_stack]{"range"}[$index]}[1];
	    }
	  else
	    {
	    $order = $print_stack[$#print_stack]{"order_hash"}{$print_stack[$#print_stack]{"array_name"}}[0];
	    $order_hash = $print_stack[$#print_stack]{"order_hash"}{$print_stack[$#print_stack]{"array_name"}}[1];
	    }
	  @$order = (grep { exists $$test_value{$_}; } (@$order));
	  if (! exists $$write_options{"limit_to_order"})
	    { push @$order, (grep { my $test = $_; ! (grep {$_ eq $test; } (@$order)); } (sort {$a cmp $b; } (keys %$test_value))); }
	  
	  push @print_stack, {"type" => "HASH", "range" => [ @$order ], "index" => 0, "original" => $test_value, "order_hash" => $order_hash }; 
	  }
	else
	  { push @print_stack, {"type" => "HASH", "range" => [ sort {$a cmp $b } (keys %$test_value) ], "index" => 0, "original" => $test_value }; }
	$hash_depth++;
	}
      elsif ($test_ref eq "ARRAY") # test value is an array
	{ 
	if (exists $print_stack[$#print_stack]{"order_hash"} && (ref($print_stack[$#print_stack]{"order_hash"}) eq "HASH"))
	  { push @print_stack, {"type" => "ARRAY", "range" => $test_value, "index" => 0, "original" => $test_value, "array_name" => $print_stack[$#print_stack]{"range"}[$index], "order_hash" => $print_stack[$#print_stack]{"order_hash"} }; }
	else
	  { push @print_stack, {"type" => "ARRAY", "range" => $test_value, "index" => 0, "original" => $test_value, "array_name" => $print_stack[$#print_stack]{"range"}[$index] }; }
	print TABLEFILE "  "x($hash_depth),$print_stack[$#print_stack]{"array_name"}, " :: ";
	}
      elsif ($print_stack[$#print_stack]{"type"} eq "HASH") # print section is a hash
	{
	print TABLEFILE "  "x($hash_depth),$print_stack[$#print_stack]{"range"}[$index]," : ","\""x(($test_value =~ /[:\[\]\,\#]/) || ($test_value =~ /\s+/) || ($test_value eq "")), $test_value , "\""x(($test_value =~ /[:\[\]\,\#]/) || ($test_value =~ /\s+/) || ($test_value eq "")),"\n";
	$print_stack[$#print_stack]{"index"}++;
	}
      else # print section is an array
	{
	if ($index)
	  { print TABLEFILE " , "; }
	print TABLEFILE "\""x(($test_value =~ /[:\[\]\,\#]/) || ($test_value =~ /\s+/) || ($test_value eq "")), $test_value , "\""x(($test_value =~ /[:\[\]\,\#]/) || ($test_value =~ /\s+/) || ($test_value eq ""));
	$print_stack[$#print_stack]{"index"}++;
	}
      }
    else # pop finished print section.
      {
      if (@print_stack > 1)
	{
	if ($print_stack[$#print_stack]{"type"} eq "HASH")
	  { 
	  print TABLEFILE "  "x($hash_depth), "]"; 
	  $hash_depth--;
	  }

	print TABLEFILE "\n";
	}
      
      pop @print_stack;
      if (@print_stack)
	{ 
	$print_stack[$#print_stack]{"index"}++; 
	if (($print_stack[$#print_stack]{"type"} eq "ARRAY") && ($print_stack[$#print_stack]{"index"} < @{$print_stack[$#print_stack]{"range"}}))
	  { print TABLEFILE "  "x($hash_depth),$print_stack[$#print_stack]{"array_name"}, " :: "; }
	}
      }
    }
  
  if ($filename ne "-")
    { close TABLEFILE; }

  return 1;
  }


# getTableValue
#   Returns value in table using list of keys and array indeces.
#
# Parameters:
#	$table - ref to table.
#	$key1, $key2, ... - list of keys and array indeces.
#
sub getTableValue
  {
  my $table = shift @_;

  return &getTableValueOrDefault($table, undef, @_);
  }

# getTableValueOrDefault
#   Returns value in table using list of keys and array indeces or returns default value.
#
# Parameters:
#	$table - ref to table.
#	$default_value - default value if keys/array indeces are not present.
#	$key1, $key2, ... - list of keys and array indeces.
#
sub getTableValueOrDefault
  {
  my $table = shift @_;
  my $default_value = shift @_;

  while (@_)
    {
    my $key = shift @_;
    if (ref($table) eq "ARRAY")
      {
      if ($key =~ /^\d+$/)
	{
	if ($key < @$table)
	  { $table = $$table[$key]; }
	else
	  { return $default_value; }
	}
      else
	{ 
	$table = $$table[0]; 
	unshift @_, $key;
	}
      }
    elsif ((ref($table) eq "HASH") && (exists $$table{$key}))
      { $table = $$table{$key}; }
    else
      { return $default_value; }
    }

  if (ref($table) eq "ARRAY")
    { return $$table[0]; }

  return $table;
  }



# testTableValue
#   Returns true if value exists in table.
#
# Parameters:
#	$table - ref to table.
#	$key1, $key2, ... - list of keys and array indeces.
#
sub testTableValue
  {
  my $table = shift @_;

  while (@_)
    {
    my $key = shift @_;
    if (ref($table) eq "ARRAY")
      {
      if ($key =~ /^\d+$/)
	{
	if ($key < @$table)
	  { $table = $$table[$key]; }
	else
	  { return 0; }
	}
      else
	{ 
	$table = $$table[0]; 
	unshift @_, $key;
	}
      }
    elsif ((ref($table) eq "HASH") && (exists $$table{$key}))
      { $table = $$table{$key}; }
    else
      { return 0; }
    }

  return 1;
  }

# numTableValue
#   Returns number of values with the given key_list in the table.
#
# Parameters:
#	$table - ref to table.
#	$key1, $key2, ... - list of keys and array indeces.
#
sub numTableValue
  {
  my $table = shift @_;

  while (@_)
    {
    my $key = shift @_;
    if (ref($table) eq "ARRAY")
      {
      if ($key =~ /^\d+$/)
	{
	if ($key < @$table)
	  { $table = $$table[$key]; }
	else
	  { return 0; }
	}
      else
	{ 
	$table = $$table[0]; 
	unshift @_, $key;
	}
      }
    elsif ((ref($table) eq "HASH") && (exists $$table{$key}))
      { $table = $$table{$key}; }
    else
      { return 0; }
    }

  if (ref($table) eq "ARRAY")
    { return scalar(@$table); }

  return 1;
  }


# testTableStructure
#   Tests a table against a template structure and returns 1 if it matches.
#
# Parameters:
#	$table - ref to table to test.
#	$structure - ref to template structure.
#
sub testTableStructure
  {
  my $table = shift @_;
  my $structure = shift @_;

  if (ref($structure) ne "HASH")
    {
    if ($structure eq "@")
      {
      if (ref($table) eq "ARRAY")
	{ return 1; }
      else
	{ return 0; }
      } 

    return $structure; 
    }

  if (ref($table) eq "ARRAY")
    {
    foreach my $item (@$table)
      { 
      if (! &testTableStructure($item,$structure))
	{ return 0; }
      }
    }
  elsif (ref($table) eq "HASH")
    {
    foreach my $key (keys %$structure)
      {
      if (! exists $$table{$key} || ! &testTableStructure($$table{$key},$$structure{$key}))
	{ return 0; }
      }
    }
  else
    { return 0; }
    
  return 1;
  }

# module must return true
return 1;
