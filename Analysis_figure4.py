import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import stats

def hamming_distance(s1, s2):
    "Return the Hamming Distance (HD) between two patterns"
    if len(s1) != len(s2):
        raise ValueError("Unequal length of patterns")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1,s2))


def analysis_of_results(Trials,mossy,maindir):

    scale_fac = 4


    N_output = 500*scale_fac
    N_input  = 100*scale_fac


    Fout_ALL = []
    similarityO_ALL = []
    similarityI_ALL = []
    errors_ALL      = []
    GC_activity = []
    O_rateALL = []
    overlap = '80'

    if mossy == 'Control':
        path=maindir+'/results/Control/'
    elif mossy == 'noMC':
        path=maindir+'/results/MC_delete/'
    elif mossy == 'noMCBC':
        path=maindir+'/results/MC_to_IN_delete/'
    elif mossy == 'noMCGC':
        path=maindir+'/results/MC_to_GC_delete/'

    Fout = []
    Fin = []
    O_rate = []

    # Number of trials-runs
    for i_trial in Trials:
        GC_active = []
        I1 = np.load(path+'input_pattern0d_'+overlap+'_'+str(i_trial)+'_'+str(1)+'.npy')
        O1 = np.load(path+'output_pattern0d_'+overlap+'_'+str(i_trial)+'_'+str(1)+'.npy')

        I2 = np.load(path+'input_pattern0d_'+overlap+'_'+str(i_trial)+'_'+str(2)+'.npy')
        O2 = np.load(path+'output_pattern0d_'+overlap+'_'+str(i_trial)+'_'+str(2)+'.npy')



        for rate_i in xrange(N_output):
            O_rate.append(O1[rate_i])
            O_rate.append(O2[rate_i])
        O_rateALL.append(O_rate)

        for n_i in xrange(N_input):
            if I1[n_i] < 1:
                I1[n_i] = 0
            else:
                I1[n_i] = 1

            if I2[n_i] < 1:
                I2[n_i] = 0
            else:
                I2[n_i] = 1

        for n_o in xrange(N_output):
            if O1[n_o] < 1:
                O1[n_o] = 0
            else:
                O1[n_o] = 1

            if O2[n_o] < 1:
                O2[n_o] = 0
            else:
                O2[n_o] = 1

        GC_active.append(len(find(O1 != 0))/float(N_output))
        GC_active.append(len(find(O2 != 0))/float(N_output))

        # mean of GC neurons that are active!!!
        s_in = len(find(I1 != 0))/float(N_input)
        s_out = max(GC_active)

        GC_activity.append((len(find(O2 != 0))/float(N_output)))

        HDin = hamming_distance(I1, I2)
        HDout = hamming_distance(O1, O2)

        finput  = float(HDin)/(2*s_in*N_input)
        foutput = float(HDout)/(2*s_out*N_output)

        if foutput < 0:
            print HDout
            print GC_active
            print

        Fout.append(foutput*100)
        Fin.append(int(finput*100))

        Fout_ALL.append(Fout)
        similarityO_ALL.append(mean(Fout))
        similarityI_ALL.append(mean(Fin))
        errors_ALL.append(stats.sem(Fout))

    return similarityO_ALL, GC_activity*100

Trials = range(1,51)
case = ['Control','noMC','noMCBC','noMCGC']
maindir = os.getcwd()
mydict1={}
mydict2={}
for icase in case:
    results1,results2 = analysis_of_results(Trials,icase,maindir)
    mydict1[icase]=results1
    mydict2[icase]=results2





fs=12 # fontsize
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
plt.hist(mydict2['Control'], bins=15, histtype='stepfilled',
         weights=np.zeros_like(mydict2['Control']) + 1. / len(mydict2['Control']), color='green', alpha=0.8, label='Control')
plt.hist(mydict2['noMC'], bins=15, histtype='stepfilled',
         weights=np.zeros_like(mydict2['noMC']) + 1. / len(mydict2['noMC']), color='purple', alpha=0.8, label='Deletion')
plt.ylabel('density', fontsize=fs+4)
plt.xlabel('Fraction of GCs recruited', fontsize=fs+4)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
plt.hist(mydict2['Control'], bins=15, histtype='stepfilled',
         weights=np.zeros_like(mydict2['Control']) + 1. / len(mydict2['Control']), color='green', alpha=0.8, label='Control')
plt.hist(mydict2['noMCGC'], bins=15, histtype='stepfilled',
         weights=np.zeros_like(mydict2['noMCGC']) + 1. / len(mydict2['noMCGC']), color='purple', alpha=0.8, label='Deletion')
plt.ylabel('density', fontsize=fs+4)
plt.xlabel('Fraction of GCs recruited', fontsize=fs+4)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
plt.hist(mydict2['Control'], bins=15, histtype='stepfilled',
         weights=np.zeros_like(mydict2['Control']) + 1. / len(mydict2['Control']), color='green', alpha=0.8, label='Control')
plt.hist(mydict2['noMCBC'], bins=15, histtype='stepfilled',
         weights=np.zeros_like(mydict2['noMCBC']) + 1. / len(mydict2['noMCBC']), color='purple', alpha=0.8, label='Deletion')
plt.ylabel('density', fontsize=fs+4)
plt.xlabel('Fraction of GCs recruited', fontsize=fs+4)


fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
plt.boxplot([mydict1['Control'],mydict1['noMC']], positions = [0.5, 1],showfliers=False)
axes.set_xticklabels(['control', 'MC deletion'], fontsize=fs+4)
axes.set_xticks([0.5, 1.0])
plt.ylim(30,80)
plt.ylabel('Output distance', fontsize=fs+4)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
plt.boxplot([mydict1['Control'],mydict1['noMCGC']], positions = [0.5, 1],showfliers=False)
axes.set_xticklabels(['control', 'MC to GC \ndeletion'], fontsize=fs+4)
axes.set_xticks([0.5, 1.0])
plt.ylim(30,80)
plt.ylabel('Output distance', fontsize=fs+4)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), sharey=True)
plt.boxplot([mydict1['Control'],mydict1['noMCBC']], positions = [0.5, 1],showfliers=False)
axes.set_xticklabels(['control', 'MC to BC \ndeletion'], fontsize=fs+4)
axes.set_xticks([0.5, 1.0])
plt.ylim(30,80)
plt.ylabel('Output distance', fontsize=fs+4)
