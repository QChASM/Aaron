package Aaron::Constants;

use strict;
use Cwd qw(cwd);
use base 'Exporter';

our @EXPORT = ();
our @EXPORT_OK = ('INFO', 'MAXSTEP', 'OTHERS');
#Code Information
use constant INFO => {
    VERSION => 1.0,
    YEAR => 2018,
    LASTUPDATE => '7/05/18',
    AUTHORS => ["Yanfei Guan", "Victoria M. Ingman", "Benjamin J. Rooks", "Steven E. Wheeler"],
};

use constant MAXSTEP => {
    TS => 5,
    INT => 4,
    TS_SINGLE => 3,
};

use constant OTHERS => {
    MAX_LAUNCH_FAILED => 5,
    MAX_CYCLE => 4,
    SLEEP_TIME => 5,
    TS_LIB=>'TS_geoms',
};

