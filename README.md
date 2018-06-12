# Aaron
AARON 1.0, An Automated Reaction Optimizer for New catalysts

This is the new object-oriented AaronTools and AARON, which will be 'publicly' released soon.  Documentation for AARON will be added coinciding with a paper describing the main features of AARON and AaronTools.

Command-line versions of the main AaronTools functionality are currently being implemented, which should allow widespread use of AaronTools even among non-Perl users.

AaronTools is an object-oriented collection of Perl modules designed to:
-Build, measure, manipulate, and compare molecular structures
-Construct input files for and parse output files from Gaussian09
-Analyze data
-Submit and monitor jobs using queuing software commonly found on high-performance computing clusters

AARON is a computational toolkit to automate the QM-based geometry optimization of transition state and intermediate structures for homogeneous catalytic reactions. It is built using AaronTools.

AARON does not implement new electronic structure methods or geometry optimization algorithms.  Instead, it is essentially an interface to Gaussian09 that automates computations.  Briefly, given a set of TS and intermediate geometries for a given catalytic reaction, AARON can automatically compute analogous structures for similar catalysts and substrates.  In this way, users can screen potential catalysts and substrates for a given reaction by predicting the activity and selectivity of a given catalyst/substrate combination.
