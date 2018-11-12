#!/usr/bin/env perl

# Move through testing directories in the appropriate order

use strict;
use warnings;

use Test::More;
use Data::Dumper;

my @tests = ( 'environment_setup',
              'job_setup' );

# run each test
my @failed;
foreach my $t (@tests) {
    eval {
        chdir($t);
        my $status = system "./$t.t 1>/dev/null 2>stderr.tmp";
        push @failed, $t if ($status);
        ok( !$status, "$t" );
        if ( $status && -f 'stderr.tmp' ) {
            open ERR, '<', 'stderr.tmp';
            while ( my $e = <ERR> ) {
                diag($e);
            }
        }
        system "rm stderr.tmp" if ( -f 'stderr.tmp' );
        chdir('..');
        1;
    } or do {
        fail("Couldn't test: $t.t");
    };
    diag($@) if $@;
}

# Summary of failed tests
diag("\nFailed tests for:") if @failed;
foreach my $f (@failed) {
    diag("    $f");
}
diag("\n");
done_testing();
