#==============================================================================
# Network of Dentate gyrus based on Chavlis et al, Hippocampus 2016
#==============================================================================

# MC delete MODEL. Mossy cells are deleted
###############################################################################


import os
from brian import *
from brian.library.ionic_currents import *
from brian.library.IF import *
import numpy as np
import time


overlap = ['80']
trial_i = [1]
Trial = trial_i[0]
trial = 1

# Initial pattern
scale_fac = 4
N_input   = 100 * scale_fac
N_granule = 500 * scale_fac
N_basket  =  25 * scale_fac
N_mossy   =  20 * scale_fac
N_hipp    =  10 * scale_fac
d_input   = 0.10 # active input density

# Active pattern of neurons
os.chdir('input_patterns/scale_'+str(scale_fac)+'/d_input_0.1')
active = list(np.load('active_pattern_'+str(Trial)+'.npy'))
inactive = [x for x in xrange(N_input) if x not in active]

reinit(states = True)
clear(erase = True, all = True)

print "\nBuilding the Network... "

# CONNECTIVITY PARAMETERS

# General parameters
E_nmda     =   0     * mV       # NMDA reversal potential
E_ampa     =   0     * mV       # AMPA reversal potential
E_gaba     = -86     * mV       # GABA reversal potential
gamma      =   0.072 * mV**-1   # Mg Concentration factor
alpha_nmda =   0.5   * ms**-1   # NMDA scale factor
alpha_ampa =   1     * ms**-1   # AMPA scale factor
alpha_gaba =   1     * ms**-1   # GABA scale factor

# EC CELLS ---> GRANULE CELLS
g_ampa_eg = 0.8066 * nS
g_nmda_eg = 1.0800 * g_ampa_eg

# EC CELLS ---> HIPP CELLS
g_ampa_eh = 0.24 * nS
g_nmda_eh = 1.15 * g_ampa_eh

# SOURCE: GRANULE CELLS
# GRANULE CELLS ---> BASKET CELLS
g_ampa_gb = 0.21 * nS
g_nmda_gb = 1.50 * g_ampa_gb

# GRANULE CELLS ---> MOSSY CELLS
g_ampa_gm = 0.50 * nS
g_nmda_gm = 1.05 * g_ampa_gm

# SOURCE: MOSSY CELLS
# MOSSY CELLS ---> GRANULE CELLS
g_ampa_mg = 0.1066 * nS
g_nmda_mg = 1.0800 * g_ampa_mg

# MOSSY CELLS ---> BASKET CELLS
g_ampa_mb = 0.35 * nS
g_nmda_mb = 1.10 * g_ampa_mb

# SOURCE: BASKET CELLS
# BASKET CELLS ---> GRANULE CELLS
g_gaba_bg = 14.0 * nS

# SOURCE: HIPP CELLS
# HIPP CELLS ---> GRANULE CELLS
g_gaba_hg = 0.12 * nS
#=======================================================================================================================


#=======================================================================================================================
# INPUT CELLS (ENTORHINAL CORTEX)
os.chdir('../../..')
from poisson_input import *
rate = 45*Hz
simtime = 1000*ms
t1 = 300 * ms
t2 =  10 * ms
spiketimes = poisson_input(active, N_input, rate, simtime, t1, t2)
Input_ec = SpikeGeneratorGroup(N_input, spiketimes)
#=======================================================================================================================

# GRANULE CELLS
# Parameters
gl_g      =   2.57   * nS  # leakage conductance
El_g      = -87.00   * mV  # reversal-resting potential
Cm_g      =   0.08   * nF  # membrane capacitance
v_th_g    = -56.00   * mV  # threshold potential
v_reset_g = -74.00   * mV  # reset potential

# AMPA/NMDA/GABA Kinetics
t_nmda_decay_g = 50.0  * ms  # NMDA decay time constant
t_nmda_rise_g  =  0.33 * ms  # NMDA rise time constant
t_ampa_decay_g =  2.5  * ms  # AMPA decay time constant
t_ampa_rise_g  =  0.1  * ms  # AMPA rise time constant
t_gaba_decay_g =  6.8  * ms  # GABA decay time constant
t_gaba_rise_g  =  0.9  * ms  # GABA rise time constant

# AMPA/NMDA/GABA Model Parameters
gamma_g      = 0.04 * mV**-1 # the steepness of Mg sensitivity of Mg unblock
Mg           = 2.0  # [mM]--mili Molar - the extracellular Magnesium concentration
eta          = 0.2 # [mM**-1] -1- mili Molar **(-1) - Magnesium sensitivity of unblock
alpha_nmda_g = 2.0  * ms**-1
alpha_ampa   = 1.0  * ms**-1
alpha_gaba   = 1.0  * ms**-1

# NOISE
g_ampa_gn = 0.008 * nS
g_nmda_gn = 0.008 * nS

# AHP patrameters
tau_ahp = 45.0 * ms
g_ahp   =  0.2 * nS

# Synaptic current equations @ SOMA
eq_soma = Equations('''
I_syn_g = I_nmda_eg + I_ampa_eg + I_nmda_mg + I_ampa_mg + I_nmda_gn + I_ampa_gn + I_gaba_bg +I_gaba_hg + I_Sahp : amp
I_nmda_eg = g_nmda_eg*(vm - E_nmda)*s_nmda_eg/(1.0 + eta*Mg*exp(-gamma*vm))                                     : amp
I_ampa_eg = g_ampa_eg*(vm - E_ampa)*s_ampa_eg                                                                   : amp
s_nmda_eg                                                                                                       : 1
s_ampa_eg                                                                                                       : 1
I_nmda_mg = g_nmda_mg*(vm - E_nmda)*s_nmda_mg/(1.0 + eta*Mg*exp(-gamma*vm))                                     : amp
I_ampa_mg = g_ampa_mg*(vm - E_ampa)*s_ampa_mg                                                                   : amp
s_nmda_mg                                                                                                       : 1
s_ampa_mg                                                                                                       : 1
I_gaba_bg = g_gaba_bg*(vm - E_gaba)*s_gaba_bg                                                                   : amp
s_gaba_bg                                                                                                       : 1
I_gaba_hg = g_gaba_hg*(vm - E_gaba)*s_gaba_hg                                                                   : amp
s_gaba_hg                                                                                                       : 1
I_nmda_gn = g_nmda_gn*(vm - E_nmda)*s_nmda_gn*(1.0 + eta*Mg*exp(-gamma*vm))                                     : amp
I_ampa_gn = g_ampa_gn*(vm - E_ampa)*s_ampa_gn                                                                   : amp
s_nmda_gn                                                                                                       : 1
s_ampa_gn                                                                                                       : 1
I_Sahp                                                                                                          : amp
dI_Sahp/dt = (g_ahp*(vm-El_g)-I_Sahp)/tau_ahp                                                                   : amp
''')

# Soma equation
granule_eqs  = MembraneEquation(Cm_g)
granule_eqs += leak_current(gl_g, El_g)
granule_eqs += IonicCurrent('I = I_syn_g : amp')
granule_eqs += eq_soma


granule = NeuronGroup(N_granule, model = granule_eqs, threshold = 'vm > v_th_g',
                     reset = 'vm = v_reset_g; I_Sahp += 0.0450*nA',
                     refractory = 20 * ms, compile = True, freeze = True)

# Initialization of membrane potential
granule.vm = El_g


#Clustering of granule cells
counter = 20
N_cl = len(granule)/counter
granule_cl = {}
for gran in xrange(N_cl):
    granule_cl[gran] = granule.subgroup(counter)
#=======================================================================================================================

#=======================================================================================================================
# BASKET CELLS

# Parameters
gl_b         =  18.054  * nS # leakage conductance
El_b         = -52      * mV # reversal-resting potential
Cm_b         =   0.1793 * nF # membrane capacitance
v_th_b       = -39      * mV # threshold potential
v_reset_b    = -45      * mV # reset potential
DeltaT_b     =   2      * mV # slope factor

# Synaptic Parameters
gamma      =   0.072 * mV**-1   # Mg Concentration factor
alpha_nmda =   0.5   * ms**-1   # NMDA scale factor
alpha_ampa =   1     * ms**-1   # AMPA scale factor
alpha_gaba =   1     * ms**-1   # GABA scale factor

#AMPA/NMDA Kinetics
t_nmda_decay_b = 130.0  * ms # NMDA decay time constant
t_nmda_rise_b  =  10.0  * ms # NMDA rise time constant
t_ampa_decay_b =   4.2  * ms # AMPA decay time constant
t_ampa_rise_b  =   1.2  * ms # AMPA rise time constant

# NOISE
g_nmda_bn       =   2.5 * nS # NMDA maximum conductance
g_ampa_bn       =   3.5 * nS # AMPA maximum conductance
t_nmda_decay_bn = 130   * ms # NMDA decay time constant
t_nmda_rise_bn  =  10   * ms # NMDA rise time constant
t_ampa_decay_bn =   4.2 * ms # AMPA decay time constant
t_ampa_rise_bn  =   1.2 * ms # AMPA rise time constant


# Synaptic current equations
eq_soma_b = Equations('''
I_syn_b = I_nmda_gb + I_ampa_gb + I_nmda_mb + I_ampa_mb + I_nmda_bn + I_ampa_bn      : amp
I_nmda_gb = g_nmda_gb*(vm - E_nmda)*s_nmda_gb*(1.0/(1 + exp(-gamma*vm)*(1.0/3.57)))  : amp
I_ampa_gb = g_ampa_gb*(vm - E_ampa)*s_ampa_gb                                        : amp
s_nmda_gb                                                                            : 1
s_ampa_gb                                                                            : 1
I_nmda_mb = g_nmda_mb*(vm - E_nmda)*s_nmda_mb*(1.0/(1 + exp(-gamma*vm)*(1.0/3.57)))  : amp
I_ampa_mb = g_ampa_mb*(vm - E_ampa)*s_ampa_mb                                        : amp
s_nmda_mb                                                                            : 1
s_ampa_mb                                                                            : 1
I_nmda_bn  = g_nmda_bn*(vm - E_nmda)*s_nmda_bn*(1.0/(1 + exp(-gamma*vm)*(1.0/3.57))) : amp
I_ampa_bn  = g_ampa_bn*(vm - E_ampa)*s_ampa_bn                                       : amp
s_nmda_bn                                                                            : 1
s_ampa_bn                                                                            : 1
''')

# Brette-Gerstner
basket_eqs = Brette_Gerstner(Cm_b, gl_b, El_b, v_th_b, DeltaT_b, tauw = 100 * ms, a = .1 * nS)
basket_eqs += IonicCurrent('I = I_syn_b : amp')
basket_eqs += eq_soma_b

basket = NeuronGroup(N_basket, model = basket_eqs, threshold = 'vm > v_th_b',
                     reset = AdaptiveReset(Vr=v_reset_b, b = 0.0205*nA),
                     refractory = 2 * ms, compile = True)

# Initialization of membrane potential
basket.vm = El_b

basket_cl = {}
for bb in xrange(N_cl):
    basket_cl[bb] = basket.subgroup(1)
#=======================================================================================================================


#=======================================================================================================================
# MOSSY CELLS

# Parameters
gl_m           =   4.53   * nS      # leakage conductance
El_m           = -64      * mV      # reversal-resting potential
Cm_m           =   0.2521 * nfarad  # membrane capacitance
v_th_m         = -42      * mV      # threshold potential
v_reset_m      = -49      * mV      # reset potential
DeltaT_m       =   2      * mV      # slope factor

# Synaptic Parameters
gamma      =   0.072 * mV**-1   # Mg Concentration factor
alpha_nmda =   0.5   * ms**-1   # NMDA scale factor
alpha_ampa =   1     * ms**-1   # AMPA scale factor
alpha_gaba =   1     * ms**-1   # GABA scale factor

#AMPA/NMDA Kinetics
t_nmda_decay_m = 100     * ms  # NMDA decay time constant
t_nmda_rise_m  =   4     * ms  # NMDA rise time constant
t_ampa_decay_m =   6.2   * ms  # AMPA decay time constant
t_ampa_rise_m  =   0.5   * ms  # AMPA rise time constant

# Noise model Parameters
g_nmda_mn       =   4.465 * nS  # NMDA maximum conductance
g_ampa_mn       =   4.7   * nS  # AMPA maximum conductance
t_nmda_decay_mn = 100     * ms  # NMDA decay time constant
t_nmda_rise_mn  =   4     * ms  # NMDA rise time constant
t_ampa_decay_mn =   6.2   * ms  # AMPA decay time constant
t_ampa_rise_mn  =   0.5   * ms  # AMPA rise time constant

# Synaptic current equations
eq_soma_m = Equations('''
I_syn_m   = I_ampa_gm + I_nmda_gm + I_ampa_mn + I_nmda_mn                           : amp
I_nmda_gm = g_nmda_gm*(vm - E_nmda)*s_nmda_gm*(1.0/(1 + exp(-gamma*vm)*(1.0/3.57))) : amp
I_ampa_gm = g_ampa_gm*(vm - E_ampa)*s_ampa_gm                                       : amp
s_nmda_gm                                                                           : 1
s_ampa_gm                                                                           : 1
I_nmda_mn = g_nmda_mn*(vm - E_nmda)*s_nmda_mn*(1.0/(1 + exp(-gamma*vm)*(1.0/3.57))) : amp
I_ampa_mn = g_ampa_mn*(vm - E_ampa)*s_ampa_mn                                       : amp
s_nmda_mn                                                                           : 1
s_ampa_mn                                                                           : 1
''')

# Brette-Gerstner
mossy_eqs  = Brette_Gerstner(Cm_m, gl_m, El_m, v_th_m, DeltaT_m, tauw = 180 * ms, a = 1 * nS)
mossy_eqs += IonicCurrent('I = I_syn_m : amp')
mossy_eqs += eq_soma_m

mossy = NeuronGroup(N_mossy, model = mossy_eqs, threshold = 'vm > v_th_m',
                     reset = AdaptiveReset(Vr=v_reset_m, b = 0.0829*nA),
                     refractory = 2 * ms, compile = True)

# Initialization of membrane potential
mossy.vm = El_m
#=======================================================================================================================

#=======================================================================================================================
# HIPP CELLS
# Parameters
gl_h      =   1.930  * nS # leakage conductance
El_h      = -59      * mV # reversal-resting potential
Cm_h      =  0.0584  * nF # membrane capacitance
v_th_h    = -50      * mV # threshold potential
v_reset_h = -56      * mV # reset potential
DeltaT_h  =   2      * mV # slope factor

# Synaptic Parameters
gamma      =   0.072 * mV**-1   # Mg Concentration factor
alpha_nmda =   0.5   * ms**-1   # NMDA scale factor
alpha_ampa =   1     * ms**-1   # AMPA scale factor
alpha_gaba =   1     * ms**-1   # GABA scale factor

#AMPA/NMDA Kinetics
t_nmda_decay_h = 110    * ms  # NMDA decay time constant
t_nmda_rise_h  =   4.8  * ms  # NMDA rise time constant
t_ampa_decay_h =  11.0  * ms  # AMPA decay time constant
t_ampa_rise_h  =   2.0  * ms  # AMPA rise time constant
# NOISE
g_nmda_hn       =   0.2 * nS  # NMDA maximum conductance
g_ampa_hn       =   0.2 * nS  # AMPA maximum conductance
t_nmda_decay_hn = 100   * ms  # NMDA decay time constant
t_nmda_rise_hn  =  5.0  * ms  # NMDA rise time constant
t_ampa_decay_hn = 11.0  * ms  # AMPA decay time constant
t_ampa_rise_hn  =  2.0  * ms  # AMPA rise time constant

# Synaptic current equations
eq_soma_h = Equations('''
I_syn_h   = I_nmda_eh + I_ampa_eh + I_nmda_hn + I_ampa_hn                           : amp
I_nmda_eh = g_nmda_eh*(vm - E_nmda)*s_nmda_eh*1./(1 + exp(-gamma*vm)/3.57)          : amp
I_ampa_eh = g_ampa_eh*(vm - E_ampa)*s_ampa_eh                                       : amp
s_nmda_eh                                                                           : 1
s_ampa_eh                                                                           : 1
I_nmda_hn = g_nmda_hn*(vm - E_nmda)*s_nmda_hn*(1.0/(1 + exp(-gamma*vm)*(1.0/3.57))) : amp
I_ampa_hn = g_ampa_hn*(vm - E_ampa)*s_ampa_hn                                       : amp
s_nmda_hn                                                                           : 1
s_ampa_hn                                                                           : 1
''')

# Brette-Gerstner
hipp_eqs  = Brette_Gerstner(Cm_h, gl_h, El_h, v_th_h, DeltaT_h, tauw = 93 * ms, a = .82 * nS)
hipp_eqs += IonicCurrent('I = I_syn_h : amp')
hipp_eqs += eq_soma_h


hipp = NeuronGroup(N_hipp, model = hipp_eqs, threshold = EmpiricalThreshold(threshold = v_th_h,refractory = 3*ms),
                     reset = AdaptiveReset(Vr=v_reset_h, b = 0.015*nA), compile = True, freeze = True)

# Initialization of membrane potential
hipp.vm = El_h
#=======================================================================================================================

#=======================================================================================================================
# ***************************************  C  O  N  N  E  C  T  I  O  N  S  ********************************************
#=======================================================================================================================
os.chdir('ConnectivityMatrices')
os.chdir('scale_'+str(scale_fac))
#  EC CELLS ----> GRANULE CELLS
a = 3.5
# Synapses at 1st dendrite
nmda_eqs_eg = '''
dj_eg/dt = -j_eg / t_nmda_decay_g + alpha_nmda_g * x_eg * (1 - j_eg) : 1
dx_eg/dt = -x_eg / t_nmda_rise_g                                     : 1
wNMDA_eg                                                             : 1
'''
synNMDA_eg = Synapses(Input_ec, granule, model = nmda_eqs_eg, pre = 'x_eg += wNMDA_eg', implicit=True, freeze=True)
granule.s_nmda_eg = synNMDA_eg.j_eg
synNMDA_eg.load_connectivity('syn_eg.txt')
synNMDA_eg.wNMDA_eg[:, :] = 1.0 * a
synNMDA_eg.delay[:, :]    = 3 * ms

ampa_eqs_eg = '''
dy_eg/dt = -y_eg / t_ampa_decay_g + alpha_ampa * h_eg * (1 - y_eg) : 1
dh_eg/dt = -h_eg / t_ampa_rise_g                                   : 1
wAMPA_eg                                                           : 1
'''
synAMPA_eg = Synapses(Input_ec, granule, model = ampa_eqs_eg, pre = 'h_eg += wAMPA_eg', implicit=True, freeze=True)
granule.s_ampa_eg = synAMPA_eg.y_eg
synAMPA_eg.load_connectivity('syn_eg.txt')
synAMPA_eg.wAMPA_eg[:, :] = 1.0 * a
synAMPA_eg.delay[:, :]    = 3 * ms


# EC CELLS ---> HIPP CELLS
# The NMDA/AMPA synapses @ hipp cell
nmda_eqs_eh = '''
dj_eh/dt = -j_eh / t_nmda_decay_h + alpha_nmda * x_eh * (1 - j_eh) : 1
dx_eh/dt = -x_eh / t_nmda_rise_h                                   : 1
wNMDA_eh                                                           : 1
'''
synNMDA_eh = Synapses(Input_ec, hipp, model = nmda_eqs_eh, pre = 'x_eh += wNMDA_eh', implicit=True, freeze=True)
hipp.s_nmda_eh = synNMDA_eh.j_eh
synNMDA_eh.load_connectivity('syn_eh.txt')
synNMDA_eh.wNMDA_eh[:, :] = 1.0
synNMDA_eh.delay[:, :]    = 3.0 * ms

ampa_eqs_eh = '''
dy_eh/dt = -y_eh / t_ampa_decay_h + h_eh*alpha_ampa*(1 - y_eh) : 1
dh_eh/dt = -h_eh / t_ampa_rise_h                               : 1
wAMPA_eh                                                       : 1
'''
synAMPA_eh = Synapses(Input_ec, hipp, model = ampa_eqs_eh, pre = 'h_eh += wAMPA_eh', implicit=True, freeze=True)
hipp.s_ampa_eh = synAMPA_eh.y_eh
synAMPA_eh.load_connectivity('syn_eh.txt')
synAMPA_eh.wAMPA_eh[:, :] = 1.0
synAMPA_eh.delay[:, :]    = 3.0 * ms

#    # GRANULE CELLS ---> MOSSY CELLS
#    # The NMDA/AMPA synapses @ mossy cell
#    nmda_eqs_gm = '''
#    dj_gm/dt = -j_gm / t_nmda_decay_m + alpha_nmda * x_gm * (1 - j_gm) : 1
#    dx_gm/dt = -x_gm / t_nmda_rise_m                                   : 1
#    wNMDA_gm                                                           : 1
#    '''
#    synNMDA_gm = Synapses(granule, mossy, model = nmda_eqs_gm, pre = 'x_gm += wNMDA_gm', implicit=True, freeze=True)
#    mossy.s_nmda_gm = synNMDA_gm.j_gm
#    synNMDA_gm.load_connectivity('syn_gm.txt')
#    synNMDA_gm.wNMDA_gm[:, :] = 1.0
#    synNMDA_gm.delay[:, :]    = 1.5 * ms
#
#    ampa_eqs_gm = '''
#    dy_gm/dt = -y_gm / t_ampa_decay_m + h_gm*alpha_ampa*(1 - y_gm) : 1
#    dh_gm/dt = -h_gm / t_ampa_rise_m                               : 1
#    wAMPA_gm                                                       : 1
#    '''
#    synAMPA_gm = Synapses(granule, mossy, model = ampa_eqs_gm, pre = 'h_gm += wAMPA_gm', implicit=True, freeze=True)
#    mossy.s_ampa_gm = synAMPA_gm.y_gm
#    synAMPA_gm.load_connectivity('syn_gm.txt')
#    synAMPA_gm.wAMPA_gm[:, :] = 1.0
#    synAMPA_gm.delay[:, :]    = 1.5 * ms

# GRANULE CELLS ---> BASKET CELLS
# The NMDA/AMPA synapses @ basket cell
synNMDA_gb = {}
synAMPA_gb = {}
for gtob in xrange(N_cl):
    nmda_eqs_gb = '''
    dj_gb/dt = -j_gb / t_nmda_decay_b + alpha_nmda * x_gb * (1 - j_gb) : 1
    dx_gb/dt = -x_gb / t_nmda_rise_b                                   : 1
    wNMDA_gb                                                           : 1
    '''
    synNMDA_gb[gtob] = Synapses(granule_cl[gtob], basket_cl[gtob], model = nmda_eqs_gb, pre = 'x_gb += wNMDA_gb', implicit=True, freeze=True)
    basket_cl[gtob].s_nmda_gb = synNMDA_gb[gtob].j_gb
    synNMDA_gb[gtob].connect_random(granule_cl[gtob], basket_cl[gtob], sparseness = 1.0)
    synNMDA_gb[gtob].wNMDA_gb[:, :] = 1.0
    synNMDA_gb[gtob].delay[:, :]    = 0.8 * ms

    ampa_eqs_gb = '''
    dy_gb/dt = -y_gb / t_ampa_decay_b + h_gb*alpha_ampa*(1 - y_gb) : 1
    dh_gb/dt = -h_gb / t_ampa_rise_b                               : 1
    wAMPA_gb                                                       : 1
    '''
    synAMPA_gb[gtob] = Synapses(granule_cl[gtob], basket_cl[gtob], model = ampa_eqs_gb, pre = 'h_gb += wAMPA_gb', implicit=True, freeze=True)
    basket_cl[gtob].s_ampa_gb = synAMPA_gb[gtob].y_gb
    synAMPA_gb[gtob].connect_random(granule_cl[gtob], basket_cl[gtob], sparseness = 1.0)
    synAMPA_gb[gtob].wAMPA_gb[:, :] = 1.0
    synAMPA_gb[gtob].delay[:, :]    = 0.8 * ms


#    # MOSSY CELLS ---> GRANULE CELLS
#    # The NMDA/AMPA synapses @ granule proximal dendrite (dendrite 2)
#    # 1st branch
#    nmda_eqs_mg = '''
#    dj_mg/dt = -j_mg / t_nmda_decay_g + alpha_nmda_g * x_mg * (1 - j_mg) : 1
#    dx_mg/dt = -x_mg / t_nmda_rise_g                                     : 1
#    wNMDA_mg                                                             : 1
#    '''
#    synNMDA_mg = Synapses(mossy, granule, model = nmda_eqs_mg, pre = 'x_mg += wNMDA_mg', implicit=True, freeze=True)
#    granule.s_nmda_mg = synNMDA_mg.j_mg
#    synNMDA_mg.load_connectivity('syn_mg.txt')
#    synNMDA_mg.wNMDA_mg[:, :] = 1.0
#    synNMDA_mg.delay[:, :]    = 3.0 * ms
#
#    ampa_eqs_mg = '''
#    dy_mg/dt = -y_mg / t_ampa_decay_g + h_mg * alpha_ampa * (1 - y_mg) : 1
#    dh_mg/dt = -h_mg / t_ampa_rise_g                                   : 1
#    wAMPA_mg                                                           : 1
#    '''
#    synAMPA_mg = Synapses(mossy, granule, model = ampa_eqs_mg, pre = 'h_mg += wAMPA_mg', implicit=True, freeze=True)
#    granule.s_ampa_mg = synAMPA_mg.y_mg
#    synAMPA_mg.load_connectivity('syn_mg.txt')
#    synAMPA_mg.wAMPA_mg[:, :] = 1.0
#    synAMPA_mg.delay[:, :]    = 3.0 * ms
#
#
#
#    # MOSSY CELL ---> BASKET CELLS
#    # The NMDA/AMPA synapses @ basket cell
#    nmda_eqs_mb = '''
#    dj_mb/dt = -j_mb / t_nmda_decay_b + alpha_nmda * x_mb * (1 - j_mb) : 1
#    dx_mb/dt = -x_mb / t_nmda_rise_b                                   : 1
#    wNMDA_mb                                                           : 1
#    '''
#    synNMDA_mb = Synapses(mossy, basket, model = nmda_eqs_mb, pre = 'x_mb += wNMDA_mb', implicit=True, freeze=True)
#    basket.s_nmda_mb = synNMDA_mb.j_mb
#    synNMDA_mb.connect_random(mossy, basket, sparseness = 1.0)
#    synNMDA_mb.wNMDA_mb[:, :] = 1.0
#    synNMDA_mb.delay[:, :]    = 3.0 * ms
#
#    ampa_eqs_mb = '''
#    dy_mb/dt = -y_mb / t_ampa_decay_b + h_mb*alpha_ampa*(1 - y_mb) : 1
#    dh_mb/dt = -h_mb / t_ampa_rise_b                               : 1
#    wAMPA_mb                                                       : 1
#    '''
#    synAMPA_mb = Synapses(mossy, basket, model = ampa_eqs_mb, pre = 'h_mb += wAMPA_mb', implicit=True, freeze=True)
#    basket.s_ampa_mb = synAMPA_mb.y_mb
#    synAMPA_mb.connect_random(mossy, basket, sparseness = 1.0)
#    synAMPA_mb.wAMPA_mb[:, :] = 1.0
#    synAMPA_mb.delay[:, :]    = 3.0 * ms

# BASKET CELLS ----> GRANULE CELLS (INHIBITION @ soma)
# Synapses @ granule cell (soma)
syn_bg = {}
for btog in xrange(N_cl):
    gaba_eqs_bg = '''
    dz_bg/dt = -z_bg / t_gaba_decay_g + r_bg*alpha_gaba*(1 - z_bg) : 1
    dr_bg/dt = -r_bg / t_gaba_rise_g                               : 1
    w_bg                                                           : 1
    '''
    syn_bg[btog] = Synapses(basket_cl[btog], granule_cl[btog], model = gaba_eqs_bg, pre = 'r_bg += w_bg', implicit=True, freeze=True)
    granule_cl[btog].s_gaba_bg = syn_bg[btog].z_bg
    syn_bg[btog].connect_random(basket_cl[btog], granule_cl[btog], sparseness = 1.0)
    syn_bg[btog].w_bg[:, :]  = 1.0
    syn_bg[btog].delay[:, :] = 0.85 * ms

# HIPP CELLS ----> GRANULE CELLS (INHIBITION @ distal dendrite)

# Synapses at granule cell distal dendrite (0)
# Synapses @ 1st branch
gaba_eqs_hg = '''
dz_hg/dt = -z_hg / t_gaba_decay_g + alpha_gaba * r_hg * (1 - z_hg) : 1
dr_hg/dt = -r_hg / t_gaba_rise_g                                   : 1
w_hg                                                               : 1
'''
syn_hg = Synapses(hipp, granule, model = gaba_eqs_hg, pre = 'r_hg += w_hg', implicit=True, freeze=True)
granule.s_gaba_hg = syn_hg.z_hg
syn_hg.load_connectivity('syn_hg.txt')
syn_hg.w_hg[:, :]  = 1.0
syn_hg.delay[:, :] = 1.6 * ms



############################################# N O I S E ################################################################
# GRANULE CELLS
noise_g = PoissonGroup(40, 2.2*Hz)
# DISTAL
# Synapses at dend00
nmda_eqs_gn = '''
dj_gn/dt = -j_gn / t_nmda_decay_g + alpha_nmda_g * x_gn * (1 - j_gn) : 1
dx_gn/dt = -x_gn / t_nmda_rise_g                                     : 1
dy_gn/dt = -y_gn / t_ampa_decay_g + alpha_ampa * h_gn * (1 - y_gn)   : 1
dh_gn/dt = -h_gn / t_ampa_rise_g                                     : 1
w_gn                                                                 : 1
'''
syn_gn = Synapses(noise_g, granule, model = nmda_eqs_gn,
                pre = 'x_gn = w_gn; h_gn = w_gn', implicit=True, freeze=True)
granule.s_nmda_gn = syn_gn.j_gn
granule.s_ampa_gn = syn_gn.y_gn
syn_gn.connect_random(noise_g, granule, sparseness = 1.0)
syn_gn.w_gn[:, :]   = 1.0
syn_gn.delay[:, :]  = '10 * rand() * ms'

# BASKET CELLS
noise_b = PoissonGroup(20*N_basket, 3*Hz)
noise_b_cl = {}
for no in xrange(N_basket):
    noise_b_cl[no] = noise_b.subgroup(20)

# Synapses at basket cell (noise role)
syn_bn = {}
for cell0 in xrange(N_basket):
    nmda_eqs_bn = '''
    dj_bn/dt = -j_bn / t_nmda_decay_bn + alpha_nmda * x_bn * (1 - j_bn) : 1
    dx_bn/dt = -x_bn / t_nmda_rise_bn                                   : 1
    dy_bn/dt = -y_bn / t_ampa_decay_bn + alpha_ampa * h_bn * (1 - y_bn) : 1
    dh_bn/dt = -h_bn / t_ampa_rise_bn                                   : 1
    w_bn                                                                : 1
    '''
    syn_bn[cell0] = Synapses(noise_b_cl[cell0], basket_cl[cell0], model = nmda_eqs_bn,
                    pre = 'x_bn = w_bn; h_bn = w_bn')
    basket_cl[cell0].s_nmda_bn = syn_bn[cell0].j_bn
    basket_cl[cell0].s_ampa_bn = syn_bn[cell0].y_bn
    syn_bn[cell0].connect_random(noise_b_cl[cell0], basket_cl[cell0], sparseness = 1.0)
    syn_bn[cell0].w_bn[:, :]  = 1.0
    syn_bn[cell0].delay[:, :] = '10 * rand() * ms'

#    # MOSSY CELLS
#    noise = PoissonGroup(30*N_mossy, 3.8*Hz)
#    noise_cl = {}
#    mossy_cl = {}
#    for no in xrange(N_mossy):
#        noise_cl[no] = noise.subgroup(20)
#        mossy_cl[no] = mossy[no]
#
#    # Synapses at mossy cell (noise role)
#    syn_mn = {}
#    for kk in xrange(N_mossy):
#        nmda_eqs_mn = '''
#        dj_mn/dt = -j_mn / t_nmda_decay_mn + alpha_nmda * x_mn * (1 - j_mn) : 1
#        dx_mn/dt = -x_mn / t_nmda_rise_mn                                   : 1
#        dy_mn/dt = -y_mn / t_ampa_decay_mn + alpha_ampa * h_mn * (1 - y_mn) : 1
#        dh_mn/dt = -h_mn / t_ampa_rise_mn                                   : 1
#        w_mn                                                                : 1
#        '''
#        syn_mn[kk] = Synapses(noise_cl[kk], mossy_cl[kk], model = nmda_eqs_mn,
#                        pre = 'x_mn = w_mn; h_mn = w_mn')
#        mossy_cl[kk].s_nmda_mn = syn_mn[kk].j_mn
#        mossy_cl[kk].s_ampa_mn = syn_mn[kk].y_mn
#        syn_mn[kk].connect_random(noise_cl[kk], mossy_cl[kk], sparseness = 1.0)
#        syn_mn[kk].w_mn[:, :]  = 1.0
#        syn_mn[kk].delay[:, :] = '10 * rand() * ms'


# HIPP Cells
noise_h = PoissonGroup(20*N_hipp, 3*Hz)
noise_h_cl = {}
hipp_cl = {}
for no in xrange(N_hipp):
    noise_h_cl[no] = noise_h.subgroup(20)
    hipp_cl[no] = hipp[no]

# Synapses at hipp cell (noise role)
syn_hn = {}
for cell2 in xrange(N_hipp):
    nmda_eqs_hn = '''
    dj_hn/dt = -j_hn / t_nmda_decay_hn + alpha_nmda * x_hn * (1 - j_hn) : 1
    dx_hn/dt = -x_hn / t_nmda_rise_hn                                   : 1
    dy_hn/dt = -y_hn / t_ampa_decay_hn + alpha_ampa * h_hn * (1 - y_hn) : 1
    dh_hn/dt = -h_hn / t_ampa_rise_hn                                   : 1
    w_hn                                                                : 1
    '''
    syn_hn[cell2] = Synapses(noise_h_cl[cell2], hipp_cl[cell2], model = nmda_eqs_hn,
                    pre = 'x_hn = w_hn; h_hn = w_hn')
    hipp_cl[cell2].s_nmda_hn = syn_hn[cell2].j_hn
    hipp_cl[cell2].s_ampa_hn = syn_hn[cell2].y_hn
    syn_hn[cell2].connect_random(noise_h_cl[cell2], hipp_cl[cell2], sparseness = 1.0)
    syn_hn[cell2].w_hn[:, :]  = 1.0
    syn_hn[cell2].delay[:, :] = '10 * rand() * ms'
#=======================================================================================================================



#=======================================================================================================================
# MONITORING
I_S = SpikeMonitor(Input_ec)
G_S = SpikeMonitor(granule)

#=======================================================================================================================

#=======================================================================================================================
# ***************************************  S  I  M  U  L  A  T  I  O  N  S  ********************************************
#=======================================================================================================================
#Simulation run

start_timestamp = time.time()

run(t1+simtime+t2, report='text', report_period = 10 *second)

sim_duration = time.time() - start_timestamp
print "\nDuration of simulation: " + str(sim_duration)

os.chdir('results/MC_delete')
output_pattern = []
for spikes in xrange(N_granule):
    output_pattern.append(len(G_S[spikes]))
np.save('output_pattern0d_'+overlap[i_save]+'_'+str(trial_i[0])+'_'+str(trial), output_pattern)

input_pattern = []
for spikes_i in xrange(len(Input_ec)):
    input_pattern.append(len(I_S[spikes_i]))
np.save('input_pattern0d_'+overlap[i_save]+'_'+str(trial_i[0])+'_'+str(trial), input_pattern)
