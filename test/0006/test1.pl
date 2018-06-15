#!/usr/bin/perl -w
use strict; use warnings;

my $rmsd;

use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

use Data::Dumper;

use AaronTools::Catalysis;

my $name = 'catalysis';
my $cata = new AaronTools::Catalysis( name => $name );

print "$name\n";


