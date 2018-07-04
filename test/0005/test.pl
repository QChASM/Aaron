#!/usr/bin/perl -w
use strict; use warnings;

my $rmsd;

eval {
    use lib $ENV{'QCHASM'};
    
    use Data::Dumper;
    
    use AaronTools::Catalysis;
   
    my $cata = new AaronTools::Catalysis( name => 'catalysis' );
    
    my $original = new AaronTools::Catalysis( name=>'ref');

    $cata->substitute('ligand', Ph=>'Me');

    $rmsd = $cata->RMSD(ref_geo=>$original);
} or do {
    my $error = $@;

    die "Error found in code: $error\n";
};

unless ($rmsd < 0.2) {
    die "Substituted structure doesn't match with the reference\n";
}

print "Test past!\n";

