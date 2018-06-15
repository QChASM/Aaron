#!/usr/bin/perl -w
use strict; use warnings;

my $rmsd;

use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

use Data::Dumper;

use AaronTools::Catalysis;

my $cata = new AaronTools::Catalysis( name => 'catalysis' );


my $subs = ['Me', 'Et', 'Cl', 'tBu'];

my @cata = $cata->screen_subs('ligand', 24=>$subs);

for my $cata (@cata) {
    $cata->printXYZ();
}


