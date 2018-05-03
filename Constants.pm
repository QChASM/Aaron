package Constants;

use strict;
use Cwd qw(cwd);
use base 'Exporter';

our @EXPORT = ();
our @EXPORT_OK = ('INFO', 'EOS', 'ADA', 'METHOD', 'HIGH_METHOD', 'BASIS', 'HIGH_BASIS', 'LOW_METHOD', 
                  'PCM', 'BOLTZMANN', 'AMU_TO_KG', 'HART_TO_KCAL', 'ROOM_TEMPERATURE', 'LAUNCH_DIR',
                  'CUTOFF', 'MAXSTEP', 'MAX_LAUNCH_FAILED', 'MAX_CYCLE', 'TEMPLATE_JOB', 'SHORT_WALL', 
                  'SLEEP_TIME', 'TS_LIB', 'NAMES', 'STANDARD_PRESSURE', 'KB', 'PLANK', 'SPEED_OF_LIGHT', 
                  'GAS_CONSTANT');

our %EXPORT_TAGS = (
                     INFORMATION => [ 'INFO' ],
                     SYSTEM => [ 'EOS', 'ADA' ],
                     THEORY => [ 'METHOD', 'HIGH_METHOD', 'BASIS',
                                 'HIGH_BASIS', 'LOW_METHOD', 'PCM'],
                     PHYSICAL => [ 'BOLTZMANN', 'AMU_TO_KG', 'HART_TO_KCAL', 'ROOM_TEMPERATURE',
                                    'KB', 'PLANK', 'SPEED_OF_LIGHT', 'GAS_CONSTANT','STANDARD_PRESSURE'],
                     OTHER_USEFUL => [ 'LAUNCH_DIR', 'MAXSTEP', 'MAX_LAUNCH_FAILED', 'MAX_CYCLE',
                                       'SHORT_WALL', 'SLEEP_TIME', 'TS_LIB', 'NAMES'],
                     COMPARE => [ 'CUTOFF'],
                     JOB_FILE => ['TEMPLATE_JOB'],
                   );

#Code Information
use constant INFO => {
    VERSION => 1.0,
    YEAR => 2017,
    LASTUPDATE => '9/07/17',
    AUTHORS => ["Yanfei Guan", "Benjamin J. Rooks", "Steven E. Wheeler"],
    AARON_HOME => "/home/einsteinguan/bin/Aaron"  		#absolute path to Aaron directory
};

#System depedent Variables
use constant EOS => {
    N_PROCS => 8,
    WALL => 48,
    SHORT_PROCS => 4,
    SHORT_WALL => 2,
    QUEUE => 'PBS'
};

use constant ADA => {
    N_PROCS => 20,
    WALL => 24,
    SHORT_PROCS => 20,
    SHORT_WALL => 12,
    NO_GET_QUOTA => 1,
    QUEUE => 'LSF'
};

#Theoretical methods
use constant {
    METHOD => 'B97D',
    HIGH_METHOD => 'wB97XD',
    BASIS => 'TZV2d2p',
    HIGH_BASIS => 'TZV2d2p',
    LOW_METHOD => 'PM6',
    PCM => 'PCM',
};

#RMSD constants
use constant CUTOFF => {
    E_CUTOFF => 0.2,
    RMSD_CUTOFF => 0.15,
    D_CUTOFF => 0.35,
    CRASH_CUTOFF => 0.8,
};

#template.job constant
use constant TEMPLATE_JOB => {
    JOB_NAME => '$jobname',
    WALL_TIME => '$walltime',
    N_PROCS => '$numprocs',
    NODE_TYPE => '$nodetype',
};

#Physical constants
use constant BOLTZMANN => 0.001987204;     #kcal/mol
use constant KB => 1.380662E-23;           #J/K
use constant AMU_TO_KG => 1.66053886E-27;
use constant HART_TO_KCAL => 627.5095;
use constant ROOM_TEMPERATURE => 298.15;
use constant PLANK => 6.62606957E-34;
use constant SPEED_OF_LIGHT => 29979245800; #cm/s
use constant GAS_CONSTANT => 1.987204;      #cal/mol
use constant STANDARD_PRESSURE => 101317;

#Other useful constants
use constant LAUNCH_DIR => cwd;
use constant MAXSTEP => {
    TS => 5,
    INT => 4,
};

use constant MAX_LAUNCH_FAILED => 5;
use constant MAX_CYCLE => 4;
use constant SHORT_WALL => 2;
use constant SLEEP_TIME => 5;
use constant TS_LIB=>'TS_geoms';

#name constants
use constant NAMES => {
    LIG_OLD => 'ORIGIN',
    LIG_NONE => 'NONE',
    TS_LIB => 'TS_geoms',
    SUBSTRATE => 'substrate',
    LIGAND => 'ligand',
    CENTER => 'center',
}
