ONEX: ONline Exploration of Time Series datasets
================================================

ONEX is a system for exploring time series datasets.
The primary feature is precomputing similarity groups for each time interval throughout the
dataset. By generating representatives for each group, an analyst is able to query the
datasets for similar time series efficiently. 
ONEX is composed of distinct phases: 
Offline operations such as group construction or outlier search and
Online operations, such as finding k-similar intervals to a query string.


Building
--------
ONEX can be built on Windows or POSIX systems using gcc (mingw) and Make.
To build ONEX, run either

`make onex`
or
`make`

To clean the codebase, run
`make clean`


Usage
-----
To run the program, start the `onex` executable in a terminal.
ONEX provides a basic Command Line Interface; you're able to enter a series of
commands to load, save, delete, generate, manipulate, query, and investigate
datasets and grouping data. You can see an example run at the bottom of this file.


Commands
--------
Run 'h', '?', or 'help' to view this list:

"""
[  0] > ?
ONEX interactive command line interface.
Available commands:
load <path>          Load a dataset from the given file. New format only.
save <id> <path>     Save a loaded dataset to the given file. New format.
drop | rm <id>       Unload a loaded database.                           
random <N> <L> <rng> Generate a random NxL dataset with range rng.       

normalize <id>       Normalize the given dataset.                        

group <id>           (Re)Group the given dataset using current ST.           
ngroup <id>          (Re)Group the given dataset using current ST (naive).   
groupST <id> <ST>    (Re)Group the given dataset using the given ST.         
ngroupST <id> <ST>   (Re)Group the given dataset using the given ST (naive). 
gload <id> <path>    Load grouping data for the given dataset from file.     
gsave <id> <path>    Save grouping data for the given dataset to file.       
gdrop <id> |grm <id> Delete grouping data for the given dataset.             

ksim <dbid> <qid> <qN> <start> <end> <k> <strict>                        
                     Get the kSimilar sequences from dataset dbid to the 
                     query string in dataset qid at qN from start to end.
                     Print k results and if strict, only check the same  
                     intervals.                                          
outlier <id> <s> <e> Find a dominant outlier in the dataset.             

describe | desc <id> Describe a dataset.                                 
print <id>           Print a dataset.                                    
list | ls            List all loaded / available datasets.               

setST <value>        Set the default Similarity Threshold.               
getST                View the current default Similarity Threshold.      

debug <value>        Set the verbosity. 0=quiet, 1=verbose.              

osave <id> <path>         Save a loaded dataset to the given file. Old format.            
oload <path> <N> <L> <D>  Load a dataset from a raw file, skip first D in each line.      

help | h | ?         Print this message.                                 
exit|quit|q          Exit from ONEX CLI.                                 
Command used 0.000136 seconds.
[  1] > exit
Quitting...
"""

Attribution
-----------
The ONEX codebase includes the trillionDTW[0] methods, uses the same pruning
methods, and even uses the trillion implementation of lemire lower/upper
bound. The simple DTW method is an extended version of one found at
bytefish.de[1].


References
----------
[0] http://www.cs.ucr.edu/~eamonn/SIGKDD_trillion.pdf
[1] http://bytefish.de/blog/dynamic_time_warping/


Example run
-----------
"""
~/src/onex $ ./onex
Welcome to ONEX.
[  0] > random 10 10 10
Generating random Time Series Set with N=10, L=10, and range=10.
Dataset successfully loaded. Index: 0
Command used 0.000159 seconds.

[  1] > ls
Available data sets:
[0] UNGROUPED <random>.
Command used 5.6e-05 seconds.

[  2] > setST 1.0
Setting default ST to 1.
Command used 9.9e-05 seconds.

[  3] > desc 0
Time series set:
Name: '<random>'
Sequences: 10
Sequence length: 10
Min / max values: 1, 10
Command used 8.9e-05 seconds.

[  4] > normalize 0
Normalizing Time Series Set 0:<random>.
Command used 3.9e-05 seconds.

[  5] >  desc 0
Time series set:
Name: '<random>'
Sequences: 10
Sequence length: 10
Min / max values: 0, 1
Command used 9.7e-05 seconds.

[  6] > group 0
Grouping Time Series Set 0:<random> with ST 1.
Initializing grouping datastructure.
Grouping intervals.
Grouping base intervals.
Grouping intersected intervals.
Finished grouping with a total of 81 groups over 55 discrete intervals.
Command used 0.000435 seconds.

[  7] > ls
Available data sets:
[0]   GROUPED <random>.
Command used 3.6e-05 seconds.

[  8] > load ndata/ECG.txt
Loading Time Series Set from file 'ndata/ECG.txt'.
Dataset successfully loaded. Index: 1
Command used 0.022749 seconds.

[  9] > desc 1
Time series set:
Name: 'ndata/ECG.txt'
Sequences: 200
Sequence length: 96
Min / max values: 0, 1
Command used 5.6e-05 seconds.

[ 10] > ls
Available data sets:
[0]   GROUPED <random>.
[1] UNGROUPED ndata/ECG.txt.
Command used 4.5e-05 seconds.

[ 11] > ksim 0 1 0 7 10 1 0
Getting the kSimilar for Time Series Set 0:<random>, query string at 0 in dataset 1:ndata/ECG.txt in interval [3,10] with k=1 and strict=0.
Found most similar interval and group: 0@[3,4]
Discovered k similar points:
Series 0, interval [3,4] is at distance 0.291642.
0.555556 0.333333 
Command used 0.000296 seconds.

[ 12] > save 0 /tmp/random.txt
Saving Time Series Set 0:<random> to file '/tmp/random.txt'.
Dataset successfully saved.
Command used 0.000267 seconds.

[ 13] > gsave 0 /tmp/random.txt
Saving grouping data for Time Series Set 0:<random> to /tmp/random.txt.
Command used 0.002382 seconds.
"""