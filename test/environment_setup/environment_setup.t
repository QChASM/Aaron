#!/usr/bin/env perl

# Test environmental variables are set appropriately
# and that necessry configuration files are found

use strict;
use warnings;

use Test::More;

ok( $ENV{QCHASM},
    "QCHASM environmental variable should be set to QCHASM path." );

my $QCHASM = $ENV{QCHASM};
ok( -f "$QCHASM/Aaron/.aaronrc",
    "$QCHASM/Aaron/.aaronrc should exist for storing group-specific configuration details."
);

if ( -f "$ENV{HOME}/.aaronrc" ) {
    pass("$ENV{HOME}/.aaronrc exists");
} else {
    diag(
        "WARNING: $ENV{HOME}/.aaronrc not found. " .
		"This is optional, but useful for storing user-specific configuration details."
    );
}

done_testing();
