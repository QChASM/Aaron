#!/usr/bin/perl -w
use strict; use warnings;

my $failed_to_submit;
my $job_found;
my $failed_to_kill;
my $QCHASM = $ENV{'QCHASM'};
eval {
    use lib $ENV{'QCHASM'};
   
    use Aaron::G_Key;
    use AaronTools::JobControl qw(get_job_template submit_job findJob killJob);

    my $Gkey = new Aaron::G_Key;

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

    ($job_found) = findJob("$QCHASM/test/0007");

    killJob($job_found) if $job_found;
    sleep(5);

    $failed_to_kill = findJob("$QCHASM/test/0007");
    1
} or do {
    my $error = $@;

    die "Error found in code: $error\n";
};


if ($failed_to_submit) {
    die "Test failed. Failed to submit test job to the queue.  Check test.job for errors.\n"
}

unless ($job_found) {
    die "Test Failed. Cannot find job submitted to the queue. Please find the job and kill manually. Contact developers to debug findjob.\n";
}

if ($failed_to_kill) {
    die "Test Failed. Cannot kill a job on the queue. Please kill the job manually.\n";
}

print "Test passed!\n";
#Leave behind .job and .log file unless test passed!
system("rm -fr test.job*");
system("rm -fr test.log");



