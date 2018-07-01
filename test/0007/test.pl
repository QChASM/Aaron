#!/usr/bin/perl -w
use strict; use warnings;

my $failed_to_submit;
my $job_found;
my $failed_to_kill;
my $AARON = $ENV{'AARON'};

eval {
    use lib $ENV{'AARON'};
   
    use AaronTools::G_Key;
    use AaronTools::JobControl qw(get_job_template submit_job findJob killJob);

    my $Gkey = new AaronTools::G_Key;

    $Gkey->read_key_from_input();

    my $wall = $Gkey->{wall};
    my $n_procs = $Gkey->{n_procs};
    my $node = $Gkey->{node};

    my $template_job = get_job_template();

    $failed_to_submit = submit_job(
        com_file=>'test.com',
        walltime=>$wall,
        numprocs=>$n_procs,
        template_job=>$template_job,
        node=>$node);

    sleep(10);

    ($job_found) = findJob("$AARON/test");

    killJob($job_found) if $job_found;
    sleep(5);

    $failed_to_kill = findJob("$AARON/test");
    1
} or do {
    my $error = $@;

    die "Error found in code: $error\n";
};

system("rm -fr test.job*");
system("rm -fr test.log");

if ($failed_to_submit) {
    die "Test failed. Failed to submit test job to the queue.\n"
}

unless ($job_found) {
    die "Test Failed. Cannot find job submitted to the queue. Please find the job and kill manually first.\n";
}

if ($failed_to_kill) {
    die "Test Failed. Cannot kill the job on the queue. Please kill the job manually first.\n";
}

print "Test passed!\n";



