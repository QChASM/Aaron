#!/usr/bin/perl -w
use strict; use warnings;

my $fail;

print "Testing for required environmental variables and files...";
unless ($ENV{'AARON'}) {
    warn "FATAL! Environmental variable \$AARON is not set.\n" .
         "Please set the path to AARON directory as \$AARON environmental variable.\n";
    $fail = 1;
}

unless ($ENV{'PERL_LIB'}) {
    warn "FATAL! Environmental variable \$PERL_LIB is not set.\n" .
         "Please set the path to your perl library as \$PERL_LIB environmental variable.\n";
    $fail = 1;
}

if ($ENV{'AARON'}) {
    my $AARON = $ENV{'AARON'};

    unless (-f "$AARON/template.job") {
        warn "WARN! No template.job found in $AARON\n" .
             "Please set up a template.job following the tutorial.\n";
        $fail = 1;
     }

     unless (-f "$AARON/.aaronrc") {
        warn "WARN! No .aaronrc file found in $AARON for the group.\n" .
              "You are highly recommanded to make a custom .aaronrc file at $AARON.\n" .
              "Otherwise, test 0007 will fail.\n";
        $fail = 1;
     }
}

my $HOME = $ENV{'HOME'};
unless (-f "$HOME/.aaronrc") {
    warn "No custom keywords file .aaronrc found in your home directory.\n".
          "You are highly recommanded to make a custom .aaronrc file at $HOME for yourself.\n";
    $fail = 1;
}

if ($fail) {
    die "Key file missed!\n"
}

print "Test past!\n";
             


