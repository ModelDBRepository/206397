import os
from brian import *
from brian.library.ionic_currents import *
from brian.library.IF import *
import numpy as np
import random as pyrandom

def input_patterns(trial_i):
    reinit(states = True)
    clear(erase = True, all = True)
    os.chdir('input_patterns')

    Trial = trial_i[0]
    # Initial pattern
    scale_fac = 2
    os.chdir('scale_'+str(scale_fac))
    N_input   = 100 * scale_fac
    d_input   = 0.10 # active input density
    os.chdir('d_input_'+str(d_input))
    # Active pattern of neurons
    active   = sorted(pyrandom.sample(xrange(N_input), int(d_input*N_input)))
    np.save('active_pattern_'+str(Trial)+'.npy', active)

    return

jobidx = int(sys.argv[1])
results = input_patterns([jobidx]) # launches multiple processes
