package AaronTools::Constants;

use strict;
use Cwd qw(cwd);
use base 'Exporter';

our @EXPORT = ();
our @EXPORT_OK = ('CUTOFF', 'TEMPLATE_JOB', 'PHYSICAL', 'UNIT', 'NAMES');

#RMSD constants
use constant CUTOFF => {
    E_CUTOFF => 0.2,
    RMSD_CUTOFF => 0.15,
    D_CUTOFF => 0.35,
    CRASH_CUTOFF => 0.8,
    TS_CUTOFF => 0.3,
    TS_D => 1.95,
};

#template.job constant
use constant TEMPLATE_JOB => {
    JOB_NAME => '$jobname',
    WALL_TIME => '$walltime',
    N_PROCS => '$numprocs',
    NODE_TYPE => '$nodetype',
};

#Physical constants
use constant PHYSICAL => {
    BOLTZMANN => 0.001987204,     #kcal/mol
    KB => 1.380662E-23,          #J/K
    ROOM_TEMPERATURE => 298.15,
    PLANK => 6.62606957E-34,
    SPEED_OF_LIGHT => 29979245800, #cm/s
    GAS_CONSTANT => 1.987204,      #cal/mol
    STANDARD_PRESSURE => 101317,
};

#Unit conversion
use constant UNIT => {
    AMU_TO_KG => 1.66053886E-27,
    HART_TO_KCAL => 627.5095,
};

#name constants
use constant NAMES => {
    LIG_OLD => 'ORIGIN',
    LIG_NONE => 'NONE',
    TS_LIB => 'TS_geoms',
    SUBSTRATE => 'substrate',
    LIGAND => 'ligand',
    CENTER => 'center',
}
