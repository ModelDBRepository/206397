Author: Spyridon Chavlis, schavlis [at] imbb.forth.gr

Code for "In vivo Imaging of Dentate Gyrus Mossy Cells in Behaving Mice", Danielson et al., 2017

This code replicates the results Figure 4 of the aforementioned paper.

1. connectivityMatrices.py constructs the network and saves the connectvities in ConnectivityMatrices folder. scale_frac represents the scale of the network

from console run: python connectivityMatrices.py

2. input_patterns.py creates the inputs to the network. You should run at least 50 times (50 different trials) in order to replicate the results

from console run: for i in $(seq 1 50);do python input_patterns.py $i;done


Experiment: In line with Chavlis et al., 2017 we run an INITIAL pattern and an overlapping one 80% overlap.

a. DG_INITIAL.py and DG_80.py --> Control
b. DG_INITIAL_noMC.py and DG_80_noMC.py --> Mossy Cells full deletion
c. DG_INITIAL_noMCBC.py and DG_80_noMCBC.py --> Mossy Cell to Basket Cells deletion
d. DG_INITIAL_noMCGC.py and DG_80_noMCGC.py --> Mossy Cell to Granule Cell deletion

As the model was run under an HPC cluster you should run the experiments at least 50 times. Parameter to change from run to run is trial_i = [1],trial_i = [2] etc.


The results are stored under results folder in an approipriate location that the code creates.


Analysis of the results

python Analysis_figure4.py

