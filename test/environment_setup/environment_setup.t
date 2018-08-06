#!/usr/bin/perl -w

# Test that necessary configuration files are found

use strict;
use warnings;

use Test::More;

my $QCHASM = $ENV{QCHASM};
$QCHASM =~ s/(.*)\/?$/$1/;
ok( -f "$QCHASM/AaronTools/template.job",
    "$QCHASM/AaronTools/template.job should exist (see tutorial)." );

ok( $ENV{QUEUE_TYPE},
    "QUEUE_TYPE environmental variable should be set to the queue type" );

done_testing();

