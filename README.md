## Introduction
This is the source code of my course project "*Efficient Centrality Computation and Social Influence Analysis on Multiplex Networks*". Please refer to the paper for the algorithm details.

## Algorithms

Three algorithms are included in this implementation code.

1. CCMN: the baseline.
1. SCCMN: our basic algorithm.
1. PSCCMN: parallel algorithm of SCCMN.


## Data
We test our code on two real-world multiplex social networks: MoscowAthletics2013 and NYClimateMarch2014. These datasets can be download from [Manlio De Domenico's homepage](https://manliodedomenico.com/data.php).



## Dependency

1. `g++` and `boost`

All the codes are runnable after cmake and make: `cmake -DCMAKE_BUILD_TYPE=Release ./` and `make -j`.
