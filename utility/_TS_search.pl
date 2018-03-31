#!/usr/bin/perl -w

use strict;
use lib $ENV{'AARON'};
use lib $ENV{'PERL_LIB'};

my $AARON = $ENV{'AARON'};

use G09Job;
use AaronInit qw($template_job);
use AaronTools::Geometry;
use AaronOutput qw(init_log);

use Cwd qw(cwd);
use Constants qw(:OTHER_USEFUL :SYSTEM :JOB_FILE :PHYSICAL);
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my $system_ada = ADA;

my $parent = cwd;
my $system = { WALL => $system_ada->{WALL},
               N_PROCS => $system_ada->{N_PROCS} };

my $Gkey = { charge => 0,
             mult => 1,
             level => new AaronInit::Theory_level(),
             solvent => 'gas',
             pcm => 'pcm',
             temperature => ROOM_TEMPERATURE };

my $Ckey = { sleeptime => SLEEP_TIME,
             parent => $parent};

my @basis;
my $method;
my $ecp;
my @chargemult;

#read arguments
GetOptions(
    'restart' => \$Ckey->{restart},
    'record' => \$Ckey->{record},
    'method|m=s' => \$method,
    'wall|w=i' => \$system->{WALL},
    'process|p=i' => \$system->{N_PROCS},
    'node|n=i' => \$system->{node},
    'basis|b=s' => \@basis,
    'ecp|e=s' => \$ecp,
    'sleep=s' => \$Ckey->{sleeptime},
    'denfit' => \$Gkey->{denfit},
    'chargemult|c=i{2}' => \@chargemult,
    'solvent|s=s' => \$Gkey->{solvent},
    'pcm=s' => \$Gkey->{pcm},
);

$Gkey->{charge} = shift @chargemult if $chargemult[0];
$Gkey->{mult} = shift @chargemult if $chargemult[0];

$system->{SHORT_WALL} = $system->{WALL};
$system->{SHORT_PROCS} = $system->{N_PROCS};

my ($input_xyz) = grep { $_ =~ /\.xyz$/ } @ARGV;
my ($input_name) = $input_xyz =~ /(\S+)\.xyz/;

my $geometry = new AaronTools::Geometry( name => $input_name );

$Gkey->{level}->read_method($method);
@basis && do{$Gkey->{level}->read_basis($_) for (@basis)};
$ecp && do{$Gkey->{level}->read_ecp($ecp)};

#read gen
$Gkey->{level}->check_gen();

my $G09job = new G09Job_TS_Single(
    name => $input_name,
    catalysis => $geometry,
    Gkey => $Gkey,
    Ckey => $Ckey,
    system => $system,
    template_job => $template_job );

#main 
init_log($input_name);
$G09job->build_com( directory => '.');
$G09job->set_status('2submit');

while ($G09job->job_running()) {
    $G09job->check_status_run();
}



