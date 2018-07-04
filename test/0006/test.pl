#!/usr/bin/perl -w
use strict; use warnings;

my @rmsd;
my @subs = ('Me', 'Et', 'Cl', 'tBu');

eval {
    use lib $ENV{'QCHASM'};
    
    use Data::Dumper;
    
    use AaronTools::Catalysis;
   
    my $cata = new AaronTools::Catalysis( name => 'catalysis' );

    my @refs;

    for my $sub (@subs) {
        my $ref = new AaronTools::Catalysis(name=>$sub);
        push (@refs, $ref);
    }

    my @catas = $cata->screen_subs('ligand', 14=>[@subs]);

    for my $i (0..3) {
        my $rmsd = $catas[$i]->RMSD(ref_geo=>$refs[$i]);
        push (@rmsd, $rmsd);
    }
    1
} or do {
    my $error = $@;

    die "Error found in code: $error\n";
};

my @fails;
for my $i (0..3) {
    unless ($rmsd[$i] < 0.2) {
        push (@fails, $subs[$i]);
    }
}

if (@fails) {
    my $fail_string = join(', ', @fails);
    die "Substitutions of $fail_string don't match with the reference structure.\n";
} 

print "Test past!\n";

