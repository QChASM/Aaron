# Aaron
AARON 1.0, An Automated Reaction Optimizer for New catalysts

We are currently merging several development branches of AARON and AaronTools, to be released soon as AARON 1.0 along with accompanying documentation.

Command-line versions of the main AaronTools functionality have been implemented, which should allow widespread use of AaronTools even among non-Perl users.

AaronTools is an object-oriented collection of Perl modules designed to:

-Build, measure, manipulate, and compare molecular structures

-Construct input files for and parse output files from Gaussian09

-Analyze data

-Submit and monitor jobs using queuing software commonly found on high-performance computing clusters

AARON is a computational toolkit to automate the QM-based geometry optimization of transition state and intermediate structures for homogeneous catalytic reactions. It is built using AaronTools.

AARON does not implement new electronic structure methods or geometry optimization algorithms.  Instead, it is essentially an interface to Gaussian09 that automates computations.  Briefly, given a set of TS and intermediate geometries for a given catalytic reaction, AARON can automatically compute analogous structures for similar catalysts and substrates.  In this way, users can screen potential catalysts and substrates for a given reaction by predicting the activity and selectivity of a given catalyst/substrate combination.
