#
#Aaron::AaronTools     Tools for working with molecular structures and G09 input and output files
#
#
# Copyright (c) 2017 Steven E. Wheeler and Yanfei Guan. All rights reserved. This program
# is free software; you can redistribute it and/or modify it under the
# same terms as Perl itself.
#

package AaronTools;
#Perl module containing commonly used subroutines for building molecules, manipulating coordinates, parsing G09 log files, interacting with LSF and PBS queues, etc.

use strict;
use warnings;
use lib $ENV{'PERL_LIB'};

use Constants qw(:INFORMATION :SYSTEM :PHYSICAL :JOB_FILE);

use Cwd;
use Switch;
use Math::Trig;
use File::Basename;
use List::Util qw(min max);
use Math::Vector::Real;
use Math::MatrixReal;

use Exporter qw(import);
our @EXPORT = qw($AARON $queue_type $LSF %radii %vdw_radii %masses finished grab_coords get_constraint mirror_coords same_structure copy_coords combine_coords printXYZ printPDB write_com num_imaginary grab_freqs print_freqs find_str findJob submit make_job_file change_distance distance angle dihedral change_dihedral set_dihedral get_quota coord_shift remove_atoms remove_fragment get_error get_connected get_connectivity check_connectivity rotate genrotate quat_rot center_genrotate center_ring substitute fused_ring shortest_path get_all_connected get_gradient get_energy get_homo_lumo get_thermo RMSD_align minimize_torsion LJ_energy calculate_ee match_molecules sub_geo map_catalyst clean_structure fix_coords get_molecules sub_rotate RMSD random_string flatten read_IRC nt_builder);

#get system information from environement vironment
our $queue_type = $ENV{'QUEUE_TYPE'} ? 
                 $ENV{'QUEUE_TYPE'} : ADA->{QUEUE};						#flag for queue type
our $g09root = $ENV{'G09_ROOT'} ?
              $ENV{'G09_ROOT'} : "/software/lms/g09_D01";       #absolute path to root directory for Gaussian09
our $AARON = $ENV{'AARON'} ?
             $ENV{'AARON'} : "/home/einsteinguan/bin/AARON";

#Physical consants
my $boltzmann = 0.001987204;			#Boltzmann constant in kcal/mol
my $h = 6.62606957E-34; 			#J.s
my $kb = 1.380662E-23; 				#J/K
my $c = 29979245800; 				#cm/s
my $R = 1.987204; 				#cal/mol
my $kcal2hartree = 0.0015936;
my $amu2kg = 1.66053886E-27;
my $hart2kcal = 627.5095;
my $ang2bohr = 1.889725989;

#ATOM data
#covalent radii (from jmol source code)
our %radii = ( H=>0.32, He=>0.93, Li=>1.23, Be=>0.90, B=>0.82, C=>0.77, N=>0.75, O=>0.73, F=>0.72, Ne=>0.71, Na=>1.54, Mg=>1.36, Al=>1.18, Si=>1.11, P=>1.06, S=>1.02, Cl=>0.99, Ar=>0.98, K=>2.03, Ca=>1.74, Sc=>1.44, Ti=>1.32, V=>1.22, Cr=>1.18, Mn=>1.17, Fe=>1.17, Co=>1.16, Ni=>1.15, Cu=>1.17, Zn=>1.25, Ga=>1.26, Ge=>1.22, As=>1.20, Se=>1.16, Br=>1.14, Kr=>1.12, Rb=>2.16, Sr=>1.91, Y=>1.62, Zr=>1.45, Nb=>1.34, Mo=>1.30, Tc=>1.27, Ru=>1.25, Rh=>1.25, Pd=>1.28, Ag=>1.34, Cd=>1.48, In =>1.44, Sn=>1.41, Sb=>1.40, Te=>1.36, I=>1.33, Xe=>1.31, Cs=>2.35, Ba=>1.98, La=>1.69, Lu=>1.60, Hf=>1.44, Ta=>1.34, W=>1.30, Re=>1.28, Os=>1.26, Ir=>1.27, Pt=>1.30, Au=>1.34, Hg=>1.49, Tl=>1.48, Pb=>1.47, Bi=>1.46, X=>0);

our %vdw_radii = (H => 1.20, He => 1.40, Li => 1.82, Be => 1.3725, B => 0.795, C => 1.70, N => 1.55, O => 1.52, F => 1.47, Ne => 1.54, Na => 2.27, Mg => 1.73, Al => 1.7, Si => 2.10, P => 1.80, S => 1.80, Cl => 1.75, Ar => 1.88, K => 2.75, Ca => 2.45, Sc => 1.37, Ti => 1.37, V => 1.37, Cr => 1.37, Mn => 1.37, Fe => 1.456, Co => 0.88, Ni => 0.69, Cu => 0.72, Zn => 0.74, Ga => 1.37, Ge => 1.95, As => 1.85, Se => 1.90, Br => 1.85, Kr => 2.02, Rb => 1.58, Sr => 2.151, Y => 1.801, Zr => 1.602, Nb => 1.468, Mo => 1.526, Tc => 1.360, Ru => 1.339, Rh => 1.345, Pd => 1.376, Ag => 1.27, Cd => 1.424, In => 1.663, Sn => 2.10, Sb => 2.05, Te => 2.06, I => 1.98, Xe=>2.00, Cs=>1.84, Ba=>2.243, La=>1.877, Lu=>2.17, Hf=>1.580, Ta=>1.467, W=>1.534, Re=>1.375, Os=>1.353, Ir=>1.357, Pt=>1.75, Au=>1.66, Hg=>1.55, Tl=>1.96, Pb=>2.02, Bi=>2.15, X=>0);

our %masses = (X => '0.', H => '1.00782503207', He => '4.00260325415', Li => '7.016004548', Be => '9.012182201', B => '11.009305406', C => '12.0', N => '14.00307400478', O => '15.99491461956', F => '18.998403224', Ne => '19.99244017542', Na => '22.98976928087', Mg => '23.985041699', Al => '26.981538627', Si => '27.97692653246', P => '30.973761629', S => '31.972070999', Cl => '34.968852682', Ar => '39.96238312251', K => '38.963706679', Ca => '39.962590983', Sc => '44.955911909', Ti => '47.947946281', V => '50.943959507', Cr => '51.940507472', Mn => '54.938045141', Fe => '55.934937475', Co => '58.933195048', Ni => '57.935342907', Cu => '62.929597474', Zn => '63.929142222', Ga => '68.925573587', Ge => '73.921177767', As => '74.921596478', Se => '79.916521271', Br => '78.918337087', Kr => '85.910610729', Rb => '84.911789737', Sr => '87.905612124', Y => '88.905848295', Zr => '89.904704416', Nb => '92.906378058', Mo => '97.905408169', Tc => '98.906254747', Ru => '101.904349312', Rh => '102.905504292', Pd => '105.903485715', Ag => '106.90509682', Cd => '113.90335854', In => '114.903878484', Sn => '119.902194676', Sb => '120.903815686', Te => '129.906224399', I => '126.904472681', Xe => '131.904153457', Cs => '132.905451932', Ba => '137.905247237', La => '138.906353267', Lu => '174.940771819', Hf => '179.946549953', Ta => '180.947995763', W => '183.950931188', Re => '186.955753109', Os => '191.96148069', Ir => '192.96292643', Pt => '194.964791134', Au => '196.966568662', Hg => '201.970643011', Tl => '204.974427541', Pb => '207.976652071', Bi => '208.980398734');

#Lennard-Jones C12 and C6 parameters from Autodock (http://www.csb.yale.edu/userguides/datamanip/autodock/html/Using_AutoDock_305.a.html)
my %LJC12 = ( "CC" => '2516582.400', "CN" => '1198066.249', "CO" => '820711.722', "CS" => '2905899.052', "CH" => '29108.222', "NC" => '1198066.249', "NN" => '540675.281', "NO" => '357365.541', "NS" => '1383407.742', "NH" => '10581.989', "OC" => '820711.722', "ON" => '357365.541', "OO" => '230584.301', "OS" => '947676.268', "OH" => '6035.457', "SC" => '2905899.052', "SN" => '1383407.742', "SO" => '947676.268', "SS" => '3355443.200', "SH" => '33611.280', "HC" => '29108.222', "HN" => '10581.989', "HO" => '6035.457', "HS" => '33611.280', "HH" => '81.920', "X" => '0');

my %LJC6 = ("CC" => '1228.800000', "CN" => '861.634784', "CO" => '754.059521', "CS" => '1418.896022', "CH" => '79.857949', "NC" => '861.634784', "NN" => '588.245000', "NO" => '505.677729', "NS" => '994.930149', "NH" => '48.932922', "OC" => '754.059521', "ON" => '505.677729', "OO" => '429.496730', "OS" => '870.712934', "OH" => '39.075098', "SC" => '1418.896022', "SN" => '994.930149', "SO" => '870.712934', "SS" => '1638.400000', "SH" => '92.212017', "HC" => '79.857949', "HN" => '48.932922', "HO" => '39.075098', "HS" => '92.212017', "HH" => '2.560000', "X" => '0');

#Substituent geometry, if linear 0, nonlienar 1

my $sub_geo = { 
                Me => { conf => 0 },
                CCH => { conf => 0},
                CF3 => { conf => 0 },
                CH2OH => { conf => 1,
                           deg => 2*pi/3,
                           num => 3,
                         },
                CHO => { conf => 1,
                         deg => pi,
                         num => 2,
                       },
                Cl => { conf => 0 },
                CN => { conf => 0 },
                COCH3 => { conf => 1,
                           deg => pi,
                           num => 2,
                         },
                COOH => { conf => 1,
                          deg => pi,
                          num => 2,
                        },
                F => { conf => 0 },
                H => { conf => 0 },
                iPr => { conf => 1,
                         deg => 2*pi/3,
                         num => 3,
                       },
                NCH3_2 => { conf => 1,
                            deg => 2*pi/3,
                            num => 3,
                          },
                NH2 => { conf => 1,
                         deg => 2*pi/3,
                         num => 3,
                       },
                NHCH3 => { conf => 1,
                           deg => 2*pi/3,
                           num => 3,
                         },
                NHOH => { conf => 1,
                          deg => 2*pi/3,
                          num => 3,
                        },
                NO2 => { conf => 1,
                         deg => pi/2,
                         num => 2,
                       },
                NO => { conf => 1,
                        deg => pi,
                        num => 2,
                      },
                OCF3 => { conf => 1,
                          deg => pi,
                          num => 2,
                        },
                OH => { conf => 1,
                        deg => pi,
                        num => 2,
                      },
                O => { conf => 0},
                SCH3 => { conf => 1,
                          deg => pi,
                          num => 2,
                        },
                SH => { conf => 1,
                        deg => pi,
                        num => 2,
                      },
                SiF3 => { conf => 0},
                SiH3 => { conf => 0},
                tBu => { conf => 1,
                         deg => pi/3,
                         num => 2,
                       },
                triazole_C => { conf => 1,
                                deg => pi/2,
                                num => 4,
                              },
                triazole_N => { conf =>1,
                                deg => pi/2,
                                num => 2,
                              },
                OMe => { conf => 1,
                         deg => pi,
                         num => 2,
                       },
                Ph => { conf => 1,
                        deg => pi/2,
                        num => 2,
                      },
                H => { conf => 0 },
              };

#checks if job has finished.  Returns 1 if finished, 0 else.  Accepts name of log file
sub finished {
  my $file = $_[0];
  my $norm="Normal termination";			#message for finished optimization
  open INFILE, "<$file" or die "Can't open $file!";
  while (<INFILE>) {
    if($_ =~ /$norm/) {
      close(INFILE);
      return 1;
    }
  }
  close(INFILE);
  return 0;
}  #end sub finished


sub get_constraint {
    my ($xyz) = @_;
    open (XYZ, "<$xyz") or die "cannot open $xyz:$!\n";
    my @constraints;
    while (<XYZ>) {
        /F:(\S+)/ && do { my @temp = split(/;/, $1); 
                          while (@temp) {
                              my $bond = [ split(/-/, shift @temp) ];
                              push (@constraints, $bond);
                          }
                        };
    }
    close (XYZ);
    return (@constraints);
}




#compares RMSD of two structures (and their mirror images); returns 0 or 1 if structures are unique or different 
#runs RMSD_align several times, attempting to renumber atoms after each round!
sub same_structure {
  my ($coords1_ref, $coords2_ref, $cutoff) = @_;
  my @coords1 = @{$coords1_ref};
  my @coords2 = @{$coords2_ref};
  my $rmsd = 1e99;
  ($rmsd,@coords2) = RMSD_align(\@coords1, \@coords2);
  if($rmsd < $cutoff) {
    return 1;
  }
  foreach (0..2) {
    #attempt to renumber...
    my @temp_coords2;
    foreach my $atom (0..$#coords1) {	#for each atom in coords1, find closest atom of same element in coords2 and put into temp_compare_coords
      my $short_distance = 1e99;
      my $short_atom = 0;
      foreach my $atom2 (0..$#coords1) {
        if($coords2_ref->[$atom2][0] eq $coords1[$atom][0]) {
          my $distance = distance($atom2, $atom, $coords2_ref, \@coords1);
          if($distance < $short_distance) {
            $short_distance = $distance;
            $short_atom = $atom2;
          }
        }
      }
  #push closest atom onto @temp_compare_coords and remove from @compare_coords
      push(@temp_coords2, $coords2[$short_atom]);
    }
    @coords2 = @temp_coords2;
    ($rmsd,@coords2) = RMSD_align(\@coords1, \@coords2);
    if($rmsd < $cutoff) {
      return 1;
    }
  }
  return 0;
} #end same_structure


#Prints coords in PDB format from @coords.  If some information is missing (atom type, residue, etc), just prins placeholders
#printPDB(ref_to_coords) or
#printPDB(ref_to_coords, filename)
sub printPDB {
  my ($coords_ref, $filename) = @_;
  my @coords = @{$coords_ref};

  if($filename) {
    print "Writing PDB file to $filename\n";
    open OUTFILE, ">$filename" or die "Can't open $filename";
    foreach my $atom (0..$#coords) {
      printf OUTFILE "%-6s%5d  %-3s %3s  %4s    %8.3f%8.3f%8.3f                       %s\n", $coords[$atom][5], $atom + 1, $coords[$atom][6], $coords[$atom][7], $coords[$atom][8], $coords[$atom][2], $coords[$atom][3], $coords[$atom][4], $coords[$atom][0];
    }
    print OUTFILE "END\n";
  } else {
    foreach my $atom (0..$#coords) {
      printf "%-6s%5d  %-3s %3s          %8.3f%8.3f%8.3f                       %s\n", $coords[$atom][5], $atom + 1, $coords[$atom][6], $coords[$atom][7], $coords[$atom][2], $coords[$atom][3], $coords[$atom][4], $coords[$atom][0];
    }
    print "END\n";
  }
}


#Reads G09 log file and returns number of imaginary frequencies
sub num_imaginary {
  my ($logfile) = @_;

  open INFILE, "<$logfile" or die "Can't open $logfile";
  while (<INFILE>) {
    if($_ =~ /(\d+) imaginary frequencies \(negative Signs\)/) {
      close(INFILE);
      return $1;
    }
  }
  close(INFILE);
  return 0;
} #end num_imaginary


#Frequency data stored as an array of arrays of arrays:
#my @freqs (frequency values)
#my @vectors (array of arrays of arrays)
#Later, replace with frequency object!
#grab_freqs($filename);
#returns three references: frequencies, intensities, vectors (normal mode vectors)
sub grab_freqs {
  my ($file) = @_;

  open INFILE, "<$file" or die "Can't open $file";

  my @freqs;
  my @intensities;
  my @vectors = ();
  my $numatoms;

  my $freq_num = 0;
READFILE:
  while (<INFILE>) {
    my $line = $_;
    chomp $line;
    if($line =~ /NAtoms=\s+(\d+)/) {
      $numatoms = $1;
    }
    if($line =~ / Harmonic frequencies /) {
      foreach (0..4) {	#skip down to actual frequency data
        <INFILE>;
      }
      my $next_line = "";
      while($next_line !~ / Thermochemistry /) {
        $next_line = <INFILE>;
        chomp($next_line);
        $next_line =~ s/^\s+//;
        if($next_line =~ /Frequencies/) {
          my @array = split(/\s+/, $next_line);
          shift(@array); #discard "Frequencies"
          shift(@array); #discard "--"
          my $num_freq_this_row = $#array + 1;
          push(@freqs, @array);
          # Skip over Red. masses and Frc consts
          <INFILE>;
          <INFILE>;
          $next_line = <INFILE>;
          $next_line =~ s/^\s+//;
          @array = split(/\s+/, $next_line);
          shift(@array); #discared "IR";
          shift(@array); #discared "Inten";
          shift(@array); #discared "--";
          push(@intensities, @array);
          #Skip coordinate labels
          my $next_temp = <INFILE>;
          next READFILE unless $next_temp =~ /x\s+y\s+z/i;
          foreach my $atom (0..$numatoms) {
            $next_line = <INFILE>;
            if($next_line =~ /\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
              #three frequencies 
              $vectors[$freq_num][$1-1][0] = $2;
              $vectors[$freq_num][$1-1][1] = $3;
              $vectors[$freq_num][$1-1][2] = $4;
              $vectors[$freq_num+1][$1-1][0] = $5;
              $vectors[$freq_num+1][$1-1][1] = $6;
              $vectors[$freq_num+1][$1-1][2] = $7;
              $vectors[$freq_num+2][$1-1][0] = $8;
              $vectors[$freq_num+2][$1-1][1] = $9;
              $vectors[$freq_num+2][$1-1][2] = $10;
            } elsif($next_line =~ /\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
              #two frequencies
              $vectors[$freq_num][$1-1][0] = $2;
              $vectors[$freq_num][$1-1][1] = $3;
              $vectors[$freq_num][$1-1][2] = $4;
              $vectors[$freq_num+1][$1-1][0] = $5;
              $vectors[$freq_num+1][$1-1][1] = $6;
              $vectors[$freq_num+1][$1-1][2] = $7;
            } elsif($next_line =~ /\s+(\d+)\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)/) {
              #one frequency
              $vectors[$freq_num][$1-1][0] = $2;
              $vectors[$freq_num][$1-1][1] = $3;
              $vectors[$freq_num][$1-1][2] = $4;
            }
          }
          $freq_num += $num_freq_this_row;
        }
      }
    }
  }
  return (\@freqs, \@intensities, \@vectors);
} #End sub grab_freqs

sub print_freqs {
  my ($freq_ref, $intens_ref, $vec_ref, $numatoms) = @_;
  my @frequencies = @{$freq_ref};
  my @intensities = @{$intens_ref};
  my @vectors2 = @{$vec_ref};

  foreach my $frequency (0..$#frequencies) {
    print "Freq = $frequencies[$frequency] ($intensities[$frequency])\n";
    print "X\tY\tZ\n";
    for my $atom (0..$numatoms) {
      printf "%3d %10.2f %10.2f %10.2f\n", $atom, $vectors2[$frequency][$atom][0],$vectors2[$frequency][$atom][1],$vectors2[$frequency][$atom][2];
    }
    print "\n";
  }
} #End sub grab_freq


#Finds any frequencies involving the concerted motion of a pair of atoms
#returns refences to arrays holding frequencies and intensities
#Checks to make sure normal mode vectors are in opposite directions and oriented along bond, and that change is at least X% of distance
#find_str($atom1, $atom2, \@freqs, \@intensities, \@vectors)
sub find_str {
  my ($atom1, $atom2, $freq_ref, $intens_ref, $vec_ref, $coords_ref, $debug) = @_;
  my @freqs = @{$freq_ref};
  my @intensities = @{$intens_ref};
  my @vectors = @{$vec_ref};
  my @coords = @{$coords_ref};

  my @str_freqs;
  my @str_intens;
  my $bond_vec = V($coords[$atom2][2] - $coords[$atom1][2],
                   $coords[$atom2][3] - $coords[$atom1][3],
                   $coords[$atom2][4] - $coords[$atom1][4]);
  $bond_vec /= $bond_vec->norm();

  foreach my $frequency (0..$#freqs) {
    my $vec1 = V(@{$vectors[$frequency][$atom1]});
    my $vec2 = V(@{$vectors[$frequency][$atom2]});

    my $change_vec = $vec1 - $vec2;
    my $percent_change = $change_vec->norm()/distance($atom1, $atom2, \@coords);
    if($vec1->norm() != 0 && $vec2->norm() != 0) {
      $vec1 /= $vec1->norm();
      $vec2 /= $vec2->norm();
  
      my $dot1 = $vec1 * $vec2;
      my $dot2 = $vec1 * $bond_vec;

      if($debug) {
        printf "%10.2f: %10.2f %10.2f %10.2f ", $freqs[$frequency],$dot1,$dot2,$percent_change;
      }
      if($dot1 < -0.7 && abs($dot2) > 0.7 && $percent_change > 0.2) {
        push(@str_freqs, $freqs[$frequency]);
        push(@str_intens, $intensities[$frequency]);
        if($debug) {
          print "YES";
        }
      }
      if($debug) {
        print "\n";
      }
    }
  }

  return (\@str_freqs, \@str_intens);
} #End find_str

#Returns jobIDs of all jobs (queued or running) in current directory, returns 0 if no jobs in queue match current directory
#This could be improved by searching for $Path more carefully!
#Works for PBS (default) or LSF
sub findJob {
  my $Path = $_[0];
  chomp($Path);

  #Strip leading directories off of $Path, to deal with different ways $HOME is treated
  $Path =~ s/^\S+$ENV{USER}//;
  
  my @jobIDs;

  if($queue_type =~ /LSF/i) {				#LSF queue
    my $bjobs=`bjobs -l 2> /dev/null`;
    #Combine into one line with no whitespace
    $bjobs =~ s/\s+//g;
    $bjobs =~ s/\r|\n//g;

    #First grab all jobs
    my @jobs = ($bjobs =~ m/(Job<\d+>.+?RUNLIMIT)/g);

    #parse each job looking for $Path
    foreach my $job (@jobs) {
      if ($job =~ /Job<(\d+)>\S+CWD<.+$Path>/) {
        push(@jobIDs,$1);
      }
    }
  }elsif ($queue_type =~ /PBS/i) {				#PBS (default)
    my $qstat = `qstat -fx`;
  
    #First grab all jobs
    my @jobs = ($qstat =~ m/<Job>(.+?)<\/Job>/g);
  
    #Grab jobIDs for all jobs matching $Path
    foreach my $job (@jobs) {
    	if ($job =~ m/<Job_Id>(\d+)\S+PBS_O_WORKDIR=\S+$Path</) {
    		push(@jobIDs, $1);
    	}
    }
  }elsif ($queue_type =~ /Slurm/i) {
    my @alljobs=`squeue -o %i_%Z`;
    foreach my $job (@alljobs) {
      if($job =~ /$Path/) {
        my @array = split(/_/, $job);
        push(@jobIDs,$array[0]);
      }
    }
  }
  
  if(@jobIDs) {
  	return @jobIDs;
  }
  return;
} #end sub findJob


#call as submit(com file, catalyst, walltime, nprocs, nocheck)
#Works for LSF or PBS (default)
sub submit {
  my ($com_file, $walltime, $numprocs, $template_job, $node) = @_;
  chomp(my $jobname=`basename $com_file .com`);
  chomp(my $Path = getcwd());

  #Don't submit if job is already running in current directory!
  my $job_running = 0;
  if(&findJob($Path)) {
    print "Job already running in directory of $jobname.com.  New job NOT submitted\n";
    $job_running = 1;
  }

  unless ($job_running) { 
    #if template_job were provide, use template_job to build job file
    if ($template_job->{job}) {

        my $job_file = $template_job->{job};

        my $job_found;

        my $template_pattern = TEMPLATE_JOB;

        open JOB_TEM, "<$job_file" or die "Cannot open $job_file:$!\n";
        my $job_content = do {local $/; <JOB_TEM>};
        close (JOB_TEM);

        $job_content =~ s/\Q$template_pattern->{JOB_NAME}/$jobname/g && do {$job_found = 1;};
        $job_content =~ s/\Q$template_pattern->{WALL_TIME}/$walltime/g;
        $job_content =~ s/\Q$template_pattern->{N_PROCS}/$numprocs/g;
        $job_content =~ s/\Q$template_pattern->{NODE_TYPE}/$node/g;
        #remove the formula part
        $job_content =~ s/&formula&\n(.*\n)*&formula&\n//g;

        for my $var (sort keys %{ $template_job->{formula} }) {
            my $var_value = eval($template_job->{formula}->{$var});
            $job_content =~ s/\Q$var/$var_value/g;
        }

        if ($job_found) {
            open JOB, ">$jobname.job";
            print JOB $job_content;
            close (JOB);
        }
    }
    
    unless (-e "$jobname.job") {

      my $memory=$numprocs*120000000;
      my $mb;
      if($queue_type eq 'LSF') {
        $memory = 0.8*$numprocs*2*10**9/8;		#memory in words per job
        $mb = 2700;				#memory in MB per core + padding
      }

      open JOB, ">$jobname.job";

      if($queue_type =~ /LSF/i) {					#LSF queue
        print JOB "#BSUB -J $jobname\n";
        print JOB "#BSUB -o $jobname.job.%J\n";
        print JOB "#BSUB -L /bin/bash\n";
        print JOB "#BSUB -W $walltime:00\n";
        print JOB "#BSUB -M $mb\n";
        print JOB "#BSUB -R 'rusage[mem=$mb]'\n";
        print JOB "#BSUB -n $numprocs\n";
        print JOB "#BSUB -R 'span[ptile=$numprocs]'\n";
        print JOB "export g09root=$g09root\n";
        print JOB ". \$g09root/g09/bsd/g09.profile\n";
        print JOB "trap \"rm -r \$SCRATCH/\$LSB_JOBID\" 0 1 2 3 9 13 14 15\n";
        print JOB "mkdir \$SCRATCH/\$LSB_JOBID\n";
        print JOB "cd \$SCRATCH/\$LSB_JOBID\n";
        print JOB "echo -P- $numprocs > Default.Route\n";
        print JOB "echo -M- $memory >> Default.Route\n";
        print JOB "module purge\n";
        print JOB "env\n";
        print JOB "cp \$LS_SUBCWD/*.chk .\n";
        print JOB "g09  < \$LS_SUBCWD/$jobname.com  > \$LS_SUBCWD/$jobname.log\n";
        print JOB "cp *.chk \$LS_SUBCWD/\n";
        print JOB "exit\n";
      }elsif ($queue_type =~ /PBS/i) {					#PBS (default)
        print JOB "#!/bin/bash\n";
        print JOB "#PBS -l walltime=$walltime:00:00,mem=8gb,nodes=1:ppn=$numprocs\n\n\n";
        print JOB "export g09root=$g09root\n";
        print JOB ". \$g09root/g09/bsd/g09.profile\n\n";
        print JOB "trap \"\\rm -r \$TMPDIR/\$PBS_JOBID\" 0 1 2 3 9 13 14 15\n\n";
        print JOB "mkdir \$TMPDIR/\$PBS_JOBID\n";
        print JOB "cd \$TMPDIR/\$PBS_JOBID\n\n";
        print JOB "echo -P- $numprocs > Default.Route\n";
        print JOB "echo -M- $memory >> Default.Route\n\n";
        print JOB "module purge\n\n";
        print JOB "env\n\n";
        print JOB "cp \$PBS_O_WORKDIR/*.chk .\n";
        print JOB "g09  < \$PBS_O_WORKDIR/$jobname.com > \$PBS_O_WORKDIR/$jobname.log\n\n";
        print JOB "cp *.chk \$PBS_O_WORKDIR/\n\n";
        print JOB "exit";
      }elsif ($queue_type =~ /Slurm/i) {
        print JOB "#!/bin/bash\n";
        print JOB "#\n";
        print JOB "#SBATCH -J $jobname -e $jobname.job.e%j -o $jobname.job.o%j\n";
        if ($node) {
          print JOB "#SBATCH -p $node\n";
        }else {
          print JOB "#SBATCH -p medium\n";
        }
        print JOB "#SBATCH -t $walltime:00:00 -n $numprocs --mem=56G\n";
        print JOB "#SBATCH --ntasks-per-node=$numprocs\n";
        print JOB "\n";
        print JOB "\n";
        print JOB "cd \$TMPDIR\n";
        print JOB "\n";
        print JOB "df -h\n";
        print JOB "export g09root=/sw/group/lms/sw/g09_D01\n";
        print JOB ". $g09root/g09/bsd/g09.profile\n";
        print JOB "\n";
        print JOB "echo -P- 28 > Default.Route \n";
        print JOB "echo -M- 56GB >> Default.Route \n";
        print JOB "\n";
        print JOB "module purge \n";
        print JOB "\n";
        print JOB "env\n";
        print JOB "\n";
        print JOB "g09  <  \$SLURM_SUBMIT_DIR/$jobname.com  >  \$SLURM_SUBMIT_DIR/$jobname.log\n";
        print JOB "\n";
        print JOB "df -h\n";
        print JOB "\n";
        print JOB "ls -al\n";
        print JOB "\n";
        print JOB "\cp *.wf? *.47 \$LS_SUBCWD\n";
        print JOB "\n";
        print JOB "exit\n";
      }
      close(JOB);
    }

    #Alert user if qsub (or bsub) returns error
    #FIXME
    if($queue_type eq 'LSF') {
      if(system("bsub < $jobname.job >& /dev/null")) {
        print "Submission denied!\n";
        return 1;
      }
    } else {
      if(system("qsub $jobname.job -N $jobname >& /dev/null")) { 
        print "Submission denied!\n";
        return 1;
      }
    }
  }
  return 0;
} #end sub submit


#Makes job file for arbitrary number of com files (LSF)
sub make_job_file {
  my ($jobname, $walltime, $mb, $numprocs, $memory, @comfiles) = @_;

  open JOB, ">$jobname.job" or die "Can't open $jobname.job";
  print JOB "#BSUB -J $jobname\n";
  print JOB "#BSUB -o $jobname.job.%J\n";
  print JOB "#BSUB -L /bin/bash\n";
  print JOB "#BSUB -W $walltime:00\n";
  print JOB "#BSUB -M $mb\n";
  print JOB "#BSUB -R 'rusage[mem=$mb]'\n";
  print JOB "#BSUB -n $numprocs\n";
  print JOB "#BSUB -R 'span[ptile=$numprocs]'\n";
  print JOB "export g09root=/software/lms/g09_D01\n";
  print JOB ". \$g09root/g09/bsd/g09.profile\n";
  print JOB "trap \"rm -r \$SCRATCH/\$LSB_JOBID\" 0 1 2 3 9 13 14 15\n";
  print JOB "mkdir \$SCRATCH/\$LSB_JOBID\n";
  print JOB "cd \$SCRATCH/\$LSB_JOBID\n";
  print JOB "echo -P- $numprocs > Default.Route\n";
  print JOB "echo -M- $memory >> Default.Route\n";
  print JOB "module purge\n";
  print JOB "env\n";
  print JOB "cp \$LS_SUBCWD/*.chk .\n";
  foreach (@comfiles) {
    print JOB "g09  < \$LS_SUBCWD/$_.com  > \$LS_SUBCWD/$_.log\n";
  }
  print JOB "cp *.chk \$LS_SUBCWD/\n";
  print JOB "exit\n";

  close(JOB);
} #End make_job_file


  return rad2deg($dihedral);
}

#get quota and return summary message
sub get_quota {
  my $limit;
  my $used;
  my $quota;
  my $message;

  if($queue_type eq 'LSF') {
    $quota = `/software/tamusc/local/bin/showquota`;
  } else {
    $quota = `/usr/lpp/mmfs/bin/mmlsquota`;
  }
  if($queue_type eq 'LSF') {
    if($quota =~ /scratch\s+(\S+)\s+(\S+)/) {
      ($limit, $used) = ($2, $1);
    }
    $message = sprintf "scratch used = $used of $limit";
  } else {
    if($quota =~ /scratch\s+\S+\s+(\d+)\s+(\d+)/ ) {
      ($limit, $used) = ($2, $1);
      $used /= 1048576;
      $limit /= 1048576;
      my $percent = 100*$used/$limit;
      $message = sprintf "scratch used = %.0f of %.0f GB (%0.f%%)", $used, $limit, $percent;
      if($percent > 90) {
        $message .= " (Dangerously full!)";
      }
    }
  }

  return $message;
} #end sub get_quota


sub remove_atoms {                                #remove atoms from atom list and return the removed atoms
  my ($coords_ref, @targets) = @_;
  my @removed;
  #sort targets so that atoms can be removed from last to first
  @targets = sort {$b <=> $a} @targets;

  foreach my $target (@targets) {
   my @temp = splice(@{$coords_ref}, $target, 1);
   push @removed, @temp;
  }
  return @removed;
} #End remove_atoms


#Comb through G09 log file and look for reason for failure
#returns name of error or "UNKOWN" if can't find anything obvious
sub get_error {
  my $filename = $_[0];

  my %errors = (
    "CHK"   => "NtrErr Called from FileIO",			#delete
    "EIGEN" => "Wrong number of Negative eigenvalues", 		#opt=noeigen
    "CONV"  => "Convergence failure -- run terminated.", 	#scf=qc
    "QUOTA" => "Erroneous write", 				#check quota and alert user; REMOVE error from end of file!
    "CLASH" => "Atoms too close", 				#flag as CLASH
    "CHARGEMULT" => "The combination of multiplicity",		#die and alert user to check catalyst structure or fix reaction_data!
    "REDUND" => "Bend failed for angle"                         #Using opt=cartesian
  );

  if(-e $filename) {
    open INFILE, "$filename" or die "Can't open $filename!";
    while (<INFILE>) {
      foreach my $error (keys %errors) {
        if($_ =~ /$errors{$error}/) {
          close(INFILE);
          return $error;
        } 
      }
    }
    close(INFILE);
    return "UNKNOWN";
  } 
  return "NOFILE";
} #end sub get_error


#Returns summary of the last optimization step of an optimization
#TODO: Modify to fill progress hash with numerical values between 0 and 100 based on progress towards four convergence criteria
sub get_gradient {
  my $file = $_[0];

  my $maxforce;
  my $rmsforce;
  my $maxdisp;
  my $rmsdisp;
  my @converged; 	#array to hold "YES/NO" for each of the four criteria

  if(-e $file) {
    open INFILE, "<$file";
    while (<INFILE>) {
      if($_ =~ /Threshold  Converged/) {
        my $line = <INFILE>;
        if($line =~ /Maximum Force\s+(\d\.\d+)\s+\d+\.\d+\s+(\S+)/) {
          $maxforce = $1;
          push(@converged, $2);
        }
        $line = <INFILE>;
        if($line =~ /RMS     Force\s+(\d\.\d+)\s+\d+\.\d+\s+(\S+)/) {
          $rmsforce = $1;
          push(@converged, $2);
        }
        $line = <INFILE>;
        if($line =~ /Maximum Displacement\s+(\d\.\d+)\s+\d+\.\d+\s+(\S+)/) {
          $maxdisp = $1;
          push(@converged, $2);
        }
        $line = <INFILE>;
        if($line =~ /RMS     Displacement\s+(\d\.\d+)\s+\d+\.\d+\s+(\S+)/) {
          $rmsdisp = $1;
          push(@converged, $2);
        }
      }
    }
    close INFILE;
  } else {
    return 0;		#file not found!
  }

  if(@converged) {
    return "Max Force: $maxforce ($converged[0]), RMS Force: $rmsforce ($converged[1]), Max Disp: $maxdisp ($converged[2]), RMS Disp: $rmsdisp ($converged[3])";
  } else {
    return "No steps yet...";
  }
} #end sub get_gradient

sub get_homo_lumo {
  my ($infile) = @_;
  open INFILE, $infile or die "Can't open $infile";

  my @occ;
  my @vir;

  while (<INFILE>) {
    chomp;
    if($_ =~ /The electronic state is/) {
      @occ = [];
      @vir = [];
      my $line = '';
      while ($line !~ /Condensed to atoms/) {
        $line = <INFILE>;
        chomp($line);
        if($line =~ /Alpha  occ\. eigenvalues --\s+(.+)$/) {
          push(@occ, split(/\s+/, $1));
        } elsif($line =~ /Alpha virt\. eigenvalues --\s+(.+)$/) {
          push(@vir, split(/\s+/, $1));
        }
      }
    }
  }
  return ($occ[-1], $vir[1]);
} #end get_homo_lumo

#Returns energy from G09 log file
sub get_energy {
  my ($infile) = @_;
  open INFILE, $infile or die "Can't open $infile";
  
  my $energy = 0;
  while (<INFILE>) {					#read in data from Gaussian log file
    chomp;
    if($_ =~ /SCF Done/o) {
    	my @array= split(/\s+/o,$_);
    	$energy = $array[5];
    }
  }
  close(INFILE);
  return $energy;
} #end get_energy

#returns energy, enthalpy, G, and grimme_G (grimme_G is calculated; the rest are read from teh log file)
#Based on Gaussian Thermochemistry white paper, calculate: E + ZPVE, H (0K), H(Temp), G(Temp) (http://www.gaussian.com/g_whitepap/thermo/thermo.pdf)
sub get_thermo {
  my ($infile,$T) = @_;

  my $v0 = 100; #cutoff for quasi-RRHO (in cm-1)
  my $P = 1*101317; #1 ATM
 
  my $energy; 
  my $mass; 						#molecular mass
  my $mult; 						#Spin multiplicity
  my @rottemps; 					#rotational temperatures (Kelvin)
  my @vibtemps; 					#vibrational temperatures (Kelvin)
  my @vibfreqs; 					#vibrational frequencies (cm-1)
  
  my $sigmar; 						#Rotational symmetry number
  my $ZPVE; 						#Zero point energy read from G09 output
  my $Elec_zpve; 					#Sum of electronic and zero-point Energies read from G09 output
  my $enthalpy;  					#Sum of electronic and thermal enthalpies read from G09 output
  my $G;  						#Sum of electronic and thermal Free Energies read from G09 output
  
  open INFILE, $infile or die "Can't open $infile";
  while (<INFILE>) {					#read in data from Gaussian log file
    my $line = $_;
    if($line =~ / Harmonic frequencies/) { 		#reset everything if new frequency calculation is found
      @vibtemps = ();
      @vibfreqs = ();
    }
    if($line =~ /^ Frequencies --/) {
            chomp;
      $line =~ s/^\s+//;
      my @array= split(/\s+/o,$line);
            foreach my $freq (2..$#array) {
        if($array[$freq] > 0 ) {
          push(@vibfreqs, $array[$freq]);
          push(@vibtemps, $array[$freq]*$c*$h/$kb);
        }
      }
    }
    if($_ =~ /SCF Done/o) {
    	my @array= split(/\s+/o,$_);
    	$energy = $array[5];
    }
    if($line =~ /^ Rotational constants \(GHZ\):\s+(\S+)\s+(\S+)\s+(\S+)/) {
      @rottemps = ($1, $2, $3);
      foreach my $rot (0..$#rottemps) {
        $rottemps[$rot] *= $h*(10**9)/($kb);
      }
    }
    if($line =~ /^ Molecular mass:\s+(\S+)/) {
      $mass = $1*$amu2kg;
    }
    if($line =~ /^ Sum of electronic and zero-point Energies=\s+(\S+)/) {
      $Elec_zpve = $1;
    } 
    if($line =~ /^ Sum of electronic and thermal Enthalpies=\s+(\S+)/) {
      $enthalpy = $1;
    }
    if($line =~ /^ Sum of electronic and thermal Free Energies=\s+(\S+)/) {
      $G = $1;
    }
    if($line =~ /^ Zero-point correction=\s+(\S+)/) {
      $ZPVE = $1;
    }
    if($line =~ / Multiplicity = (\d+)/) {
      $mult = $1;
    }
    if($line =~ / Rotational symmetry number\s+(\d+)/) {
      $sigmar = $1;
    }
  } #end read G09 log file
  close(INFILE);
  unless ($enthalpy) {
    return ($energy);
  }
  #Calculate average moment of inertia for Grimme's quasi-RRHO approach
  my $Bav = ($h**2/(24*pi**2*$kb))*(1/$rottemps[0] + 1/$rottemps[1] + 1/$rottemps[2]);
  
  if ($#rottemps!=2) {
    die "Problem reading Rotational constants";
  }
  
  #Translational component of Entropy
  my $qt = (2*pi*$mass*$kb*$T/($h*$h))**(3/2)*$kb*$T/$P;
  my $St = $R*(log($qt) + 5/2);
  
  #Translation component of Energy
  my $Et = 3*$R*$T/2;
  
  #Electronic component of Entropy
  my $Se = $R*(log($mult));
  
  #Rotational component of Entropy
  my $qr = (sqrt(pi)/$sigmar)*($T**(3/2)/sqrt($rottemps[0]*$rottemps[1]*$rottemps[2]));
  my $Sr = $R*(log($qr) + 3/2);
  
  #Rotational component of Energy
  my $Er = 3*$R*$T/2;
  
  #Vibrational component of Entropy and Energy
  my $Ev = 0;
  my $Sv = 0;
  my $Sv_quasiRRHO = 0;
  
  foreach my $i (0..$#vibtemps) {
    my $Sv_temp = $vibtemps[$i]/($T*(exp($vibtemps[$i]/$T)-1)) - log(1-exp(-$vibtemps[$i]/$T));
    
    $Sv += $Sv_temp;
      $Ev += $vibtemps[$i]*(1/2 + 1/(exp($vibtemps[$i]/$T) - 1));
    
    #calculate quasi-RRHO contribution to Sv
    my $mu = $h/(8*pi**2*$vibfreqs[$i]*$c);
    my $mu_prime = $mu*$Bav/($mu + $Bav);
    my $Sr = 1/2 + log(sqrt(8*pi**3*$mu_prime*$kb*$T/$h**2));
    
    my $weight = 1/(1+($v0/$vibfreqs[$i])**4);
    
    $Sv_quasiRRHO += $weight*$Sv_temp + (1-$weight)*$Sr;
  }
  
  $Sv *= $R;
  $Ev *= $R;
  $Sv_quasiRRHO *= $R;
  
  #Grab Electronic energy from $Elec_zpve and $ZPVE
  my $E_e = $Elec_zpve - $ZPVE;
  my $Etot = $Et + $Er + $Ev;
  my $Hcorr = $Etot + $R*$T;
  my $Stot = $St + $Sr + $Sv + $Se;
  my $Stot_quasiRRHO = $St + $Sr + $Sv_quasiRRHO + $Se;
  my $Gcorr = $Hcorr - $T*$Stot;
  my $Gcorr_quasiRRHO = $Hcorr - $T*$Stot_quasiRRHO;
  $Hcorr *= $kcal2hartree/1000;
  $Gcorr_quasiRRHO *= $kcal2hartree/1000;
  
  my $Grimme_G = $E_e + $Gcorr_quasiRRHO;
  return ($E_e, $enthalpy, $G, $Grimme_G);
} #end of sub get_thermo




#Given references to arrays holding (R) and (S) absolute energies, compute ee based on Boltzmann weighted sums
#my $ee = calculate_ee(\@R_energies, \@S_energies, $temp);
sub calculate_ee {
  my ($R_ref, $S_ref, $temp) = @_;

  my $RT = $boltzmann*$temp;
  my @R_vals = @{$R_ref};
  my @S_vals = @{$S_ref};
  my $R_sum = 0;
  my $S_sum = 0;
  foreach (@R_vals) {
    $R_sum += exp(-$hart2kcal*($_ - $R_vals[0])/$RT);
  }
  foreach (@S_vals) {
    $S_sum += exp(-$hart2kcal*($_ - $R_vals[0])/$RT);
  }
  my $ee = 0;
  if($R_sum + $S_sum != 0) {
    $ee = ($R_sum - $S_sum)/($R_sum + $S_sum);
  }
  return $ee;
} #end calculate_ee

#Prints random string
#random_string(16, 'a'..'z') will be a random string of 16 letters from a to z
sub random_string { 
  join'', @_[ map{ rand @_ } 1 .. shift ] 
}

#takes reference to a @coords array and converts it to a scalar with escaped line returns (for ALEXANDER mostly)
#my $geom = flatten(\@coords);
sub flatten {
  my ($coords_ref) = @_;
  my @coords = @{$coords_ref};

  my $num_atoms = $#coords + 1;
  my $geometry = "$num_atoms\\\\n\\\\n";
  foreach my $atom (0..$#coords) {
    $geometry = $geometry . "$coords[$atom][0]  $coords[$atom][2]    $coords[$atom][3]   $coords[$atom][4]\\\\n";
  }

  return $geometry;
} #end flatten


#
#@key_atoms = [1,2] #central atoms
#@other_atom = [4 1 3]
#@other_atom = [5 3 4]
#overlay to similar molecules together by rotating the lingking bond of key fragment of molecules. 
#my @coords = match_molecules(\@coords1, \@coords2, \@keyatoms1, \@keyatoms2, \@bonds1, \@bonds2)
#the @bonds1 and @bonds2 are array of arrays
#@bonds1 = (
#          [1,3],
#          [3,4]
#          [3,4]
#          )
#@keyatoms = (
#            [1,2,3],
#            [4],
#            [5],
#            [5],
#            )
# Since we need 3 keyatoms for the central part and one key atom for each fragment,the order of the peripherial fragment should match that of bonds
#note: that the atom number and the linking order of bonds must be exatly the same

   
sub fix_coords {
  my ($coords, $subs, $specified) = @_;
  my @coords = @{$coords};
  my @subs = @{$subs};
  my @specified = @{$specified};
  my @connect = get_connected(\@coords);
  my @molecules = get_molecules(@connect);
  my @fix;
  for (@molecules) {
    my @molecule = @{$_};
    my $check;
    for my $i (@subs) { 
      if (grep { $i == $_ } @molecule) {
        $check = 1;
      }
    }
    if(!$check) {
      push(@fix, @molecule);
    } 
  }
  push(@fix, @specified);
  return @fix;
} 

                
sub get_molecules {
  my @connected = @_;
  my @molecules;
  my @molecules_count;
#  #remove avoid_atom from list of atoms connected to $start_atom
  my $start_atom = 0;
  for(;;) {
    my $i=0;
    my @positions = ($start_atom);
    #I can probably use something simpler than a hash here (since I'm not really using the values)
    my %visited = ( $start_atom => '0') ; #keys are numbers of the visited atoms, values are not used

    #loop until @positions is empty
    while(@positions) {
      my $position = shift(@positions);
      foreach my $atom (@{$connected[$position]}) {       #grab all atoms connected to current atom and add to queue (unless already visited)
        if($atom >= 0 && !exists $visited{$atom}) {
          push(@positions, $atom);
          $visited{$atom} = 0;
        }
      }
    }
    my @molecule = keys %visited;
    push (@molecules, \@molecule);
    push (@molecules_count, @molecule);
    while (grep { $i == $_ } @molecules_count) {
         $i++;
    }
    $start_atom= $i;
    $start_atom == $#connected+1 && last;
  }
  return (@molecules);
}


#reads log file with IRC and and returns references to three arrays (@irc, @energies, @reaction_path)
#@irc is an array of @coords-like 2-D arrays
sub read_IRC {
  my $filename = $_[0];
  my @irc;
  my @energies;
  my @reaction_path;

  if(-e $filename) {                                   #check to make sure file exists first
    open (INFILE, "<$filename") or die "Can't open $filename";
    my @elements=('Bq','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe');

    if($filename !~ /\S+\.log/) { 
      print "Expecting .log file!\n";
      return 0;
    }
    #grab each geometry from log file
    while (<INFILE>) {
      my $line=$_;
      if($line =~ / orientation:/) {
        my @coords = ();
        $line = <INFILE>;
        $line = <INFILE>;
        $line = <INFILE>;
        $line = <INFILE>;
        $line = <INFILE>;
        do {
          if($line =~ /^\s+\d+\s+(\S+)\s+\S*\s+(\S+)\s+(\S+)\s+(\S+)/) {
            my @atom = ($elements[$1], 0, $2, $3, $4);
            push(@coords, \@atom);
            $line = <INFILE>;
          }
        } while(!($line =~ /--/));
        push(@irc, \@coords);
      } elsif ($line =~ /Summary of reaction path following/) {
        $line = <INFILE>;
        $line = <INFILE>;
        $line = <INFILE>;
        do {
          if($line =~ / \d+\s+(\S+)\s+(\S+)/) {
            push(@energies, $1*$hart2kcal);
            push(@reaction_path, $2/$ang2bohr);
            $line = <INFILE>;
          }
        } while(!($line =~ /--/));
      }
    }
  }
  return (\@irc, \@energies, \@reaction_path);
} #end read_irc


#sets a given dihedral angle
#set_dihedral(atom1, atom2, atom3, atom4, new_angle (in degrees!), \@coords)
sub set_dihedral {
  my ($atom1, $atom2, $atom3, $atom4, $new_tau, $coords_ref) = @_;
  my @coords = @{$coords_ref};

  my $tau = dihedral($atom1, $atom2, $atom3, $atom4, \@coords);
  change_dihedral($atom2, $atom3, deg2rad($new_tau - $tau), \@coords);
}


#Builds carbon nanotubes or nanotube fragment of a given size and radius of curvature
#For complete (N,N)-CNT: @coords = nt_builder($N, $length)
#For fragment: @coords = nt_builder($width, $length, $radius)
#or
#@coords = nt_builder($width, $length, $radius, $angular_offset) where angular_offset is number of rings rotated about z-axis
sub nt_builder {
  my ($width, $length, $radius, $angular_offset) = @_;
  #Standard C-C bond length and C-H bond length
  my $CC = 1.415;
  my $CH = 1.08;
  
  my $fragment = 0;
  my $angle_tally;
  my @coords;
  if(not defined($angular_offset)) {
    $angular_offset = 0;
  }
  
  if(defined($radius)) {
    $fragment = 1;
  } else {
    if($width < 2) {
      die("Can't build nanotube smaller than (2,2)!\n");
    } else {
      $radius = newton($CC, $width);  #Quick Newton-Raphson solver to get nanotube radius numerically
    }
  }

  
  my $CC_angle = 2*asin($CC/(2*$radius));
  my $CC_halfangle = 2*asin($CC/(4*$radius));
  my $CC_side = $CC*sqrt(3.0)/2.0;
  
  my @atom = ("C", 0, $radius, 0, 0);
  my @atom_coords;
  push(@atom_coords, \@atom);
  
  #Offset rotation to place center of fragment normal to x-axis
  rotate('z', -($width/2+$width-1)*$CC_angle - $angular_offset*2*($CC_angle+$CC_halfangle), \@atom_coords);
  #slide to center tube/fragment along z-axis
  #slide($CC_side*($length-1)/2);
  coord_shift(0, 0, $CC_side*($length-1)/2, \@atom_coords);
  $angle_tally = 0;
  
  for(my $row=0; $row<$length; $row++) {
    if($row%2==1) {
      if($row!=$length-1 || $fragment==0) {
        @coords = combine_coords(\@coords, \@atom_coords);
      }
      rotate('z', $CC_angle+2*$CC_halfangle, \@atom_coords);
      $angle_tally += $CC_angle+2*$CC_halfangle;
    } else {
      rotate('z', $CC_halfangle, \@atom_coords);
      $angle_tally += $CC_halfangle;
    }
    for (my $ring=0; $ring<$width; $ring++) {
      if($row!=$length-1 || $row%2!=1 || $ring != $width-1 || $fragment==0) {
        @coords = combine_coords(\@coords, \@atom_coords);
      }
      rotate('z', $CC_angle, \@atom_coords);
      $angle_tally += $CC_angle;
      if($row%2!=1 || $ring != $width-1) {
        @coords = combine_coords(\@coords, \@atom_coords);
      }
      rotate('z', $CC_angle+2*$CC_halfangle, \@atom_coords);
      $angle_tally += $CC_angle+2*$CC_halfangle;
    }
  
    #Reset and shift
  #  slide(-$CC_side);
    coord_shift(0, 0, -$CC_side, \@atom_coords); 
    rotate('z', -$angle_tally, \@atom_coords);
    $angle_tally = 0;
  }
  
  #Cap open valences
  my $numCs = $#coords;
  foreach my $atom1 (0..$numCs) {
    my @vector;
    my $neighbors = 0;
    foreach my $atom2 (0..$numCs) {
      if($atom1 != $atom2) {
        if(distance($atom1, $atom2, \@coords) < $CC+0.1) {
          $neighbors++;
          $vector[0] += $coords[$atom2][2]-$coords[$atom1][2];
          $vector[1] += $coords[$atom2][3]-$coords[$atom1][3];
          $vector[2] += $coords[$atom2][4]-$coords[$atom1][4];
        }
      }
    }
    if($neighbors < 3) {
      my $norm = sqrt($vector[0]**2+$vector[1]**2+$vector[2]**2);
      foreach (@vector) {
        $_ *= $CH/$norm
      }
    
      my @Hatom = ("H", 0, $coords[$atom1][2] - $vector[0], $coords[$atom1][3] - $vector[1], $coords[$atom1][4] - $vector[2]);
      push(@coords, [@Hatom]);
    }
    if($neighbors < 2) {
      die "Dangling carbon $atom1!";
    }
    if($neighbors > 4) {
      die "Too many neighbors for atom $atom1 (radius too small to accommodate $width rings)";
    }
  }

  return @coords;
}
  
  
#Dirty Newton solver to get radius of closed CNTs
sub newton {
  my ($CC, $width) = @_;
  #Threshold for Newton-Raphson solver to get radius
  my $THRESH = 1E-10;

  my $lastradius = 3*$CC*$width/(2*pi); #starting guess from conventional formula
  my $old_gap = get_CNT_gap($lastradius, $CC, $width);
  my $radius = $lastradius + 0.01; #arbitrary step to get process started
  my $gap = get_CNT_gap($radius, $CC, $width);

  #Simple Newton solver using very crude finite difference derivative
  while(abs($gap) > $THRESH) {
    my $newradius = $radius - $gap*($radius - $lastradius)/($gap - $old_gap);
    $old_gap = $gap;
    $gap = get_CNT_gap($newradius, $CC, $width);
    $lastradius = $radius;
    $radius = $newradius;
  }
  return $radius;
}

#Function to be minimized to get the radius of closed CNT
sub get_CNT_gap {
  my ($guess, $CC, $width) = @_;
  my $value = asin($CC/(2*$guess)) + asin($CC/(4*$guess)) - pi/(2*$width);
  return $value;
}
1;
