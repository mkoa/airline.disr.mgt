# airline.disr.mgt

This repository contains Python source code of an OR model aimed at optimizing airline schedule recovery in disruption situations


## how to use the code

This repository is composed of 2 branches:
- the master branch corresponds to the complete model, optimizing Aircraft routing + Passenger routing
- the 'AC only' branch is a simpler version of the master branch, whose model optimizes Aircraft routing only

In both branches, the source codes is basically just of bunch of helper functions plus a model builing + solving procedure called 'solve'. It takes as input:
- An excel file containing instance data (basically a list of aircrafts, passengers and flights with a few attributes)
- A time start integer representing the begginning of the time window in HHMM format (e.g. 0850, 0910) - this should be the time at which the procedure is run
- A time end integer, same as above. This should be set as far in the future as you want. Beware of the computational load if this is more than ~6 hours after time start.
- A time resolution parameter, basically controlling discretization of time, in minutes - 10 will split time in chunks of 10 minutes


The excel file in the master repository is a basic example with very few aircrafts, passenger groups and flights. Bigger schedule can be dealt with by following the same basic file structure.
