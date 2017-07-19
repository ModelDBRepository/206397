from brian import *
import numpy as np
import matplotlib.pyplot as plt
import os
import os.path
from scipy.spatial.distance import *
from scipy import stats
import itertools

def hamming_distance(s1, s2):
    "Return the Hamming Distance (HD) between two patterns"
    if len(s1) != len(s2):
        raise ValueError("Unequal length of patterns")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1,s2))
local = 0
def analysis_pattern_separation_mossy(Trials, mossy, dend):
    scale_fac = 4


    N_output = 500*scale_fac
    N_input  = 100*scale_fac

    Fout_ALL = []
    similarityO_ALL = []
    similarityI_ALL = []
    errors_ALL      = []
    GC_activity = []
    O_rateALL = []
    overlap = ['80']

    for i_over in xrange(len(overlap)):
        if mossy == 'control':
            os.chdir('results/Control')
        elif mossy == 'loss':
            os.chdir('results/MC_delete')
        os.chdir('overlap_'+overlap[i_over])
        Fout = []
        Fin = []
        O_rate = []
        for i_trial in Trials:
            GC_active = []
            if mossy == 'control':
                I1 = np.load('input_pattern' +str(dend)+'d_'+overlap[i_over]+'_'+str(i_trial)+'_'+str(1)+'.npy')
                O1 = np.load('output_pattern'+str(dend)+'d_'+overlap[i_over]+'_'+str(i_trial)+'_'+str(1)+'.npy')

                I2 = np.load('input_pattern' +str(dend)+'d_'+overlap[i_over]+'_'+str(i_trial)+'_'+str(2)+'.npy')
                O2 = np.load('output_pattern'+str(dend)+'d_'+overlap[i_over]+'_'+str(i_trial)+'_'+str(2)+'.npy')
            elif mossy == 'loss':
                I1 = np.load('input_pattern' +str(dend)+'d_'+overlap[i_over]+'_'+str(i_trial)+'_'+str(1)+'.npy')
                O1 = np.load('output_pattern'+str(dend)+'d_'+overlap[i_over]+'_'+str(i_trial)+'_'+str(1)+'.npy')

                I2 = np.load('input_pattern' +str(dend)+'d_'+overlap[i_over]+'_'+str(i_trial)+'_'+str(2)+'.npy')
                O2 = np.load('output_pattern'+str(dend)+'d_'+overlap[i_over]+'_'+str(i_trial)+'_'+str(2)+'.npy')

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

            if i_over == 0:
                GC_activity.append((len(find(O1 != 0))/float(N_output)))
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

    return similarityI_ALL, similarityO_ALL, errors_ALL, GC_activity, Fout_ALL, O_rateALL

Trials = range(1, 51)

#dend = [3, 6, 12]
dend = [0]

Ndend = dend[0]
output_MOSSY  = analysis_pattern_separation_mossy(Trials, 'control', Ndend)
output_NO_MOSSY  = analysis_pattern_separation_mossy(Trials, 'loss', Ndend)


FS = 16

N = len(output_MOSSY[0])
ind = np.arange(N)  # the x locations for the groups
width = 0.25       # the width of the bars

fig1, ax = plt.subplots()

Means_MOSSY = output_MOSSY[1]
Sem_MOSSY   = output_MOSSY[2]
rects1 = ax.bar(ind, Means_MOSSY, width, color='green', yerr=Sem_MOSSY, ecolor='k')

Means_NO_MOSSY = output_NO_MOSSY[1]
Sem_NO_MOSSY   = output_NO_MOSSY[2]
rects2 = ax.bar(ind+width, Means_NO_MOSSY, width, color='lightgreen', yerr=Sem_NO_MOSSY, ecolor='k')


# add some text for labels, title and axes ticks
ax.set_ylabel('Output Similarity', fontsize=FS)
ax.set_xlabel('Input Similarity', fontsize=FS)
#ax.set_yticks(range(0, 81, 20))
ax.set_xticks(ind + width)
ax.set_xticklabels( ('80') )
#    plt.xticks(fontsize = FS-2)
#    plt.yticks(fontsize = FS-2)
ax.legend( (rects1[0], rects2[0]), ('Control', 'Mossy loss'), loc = 'up right' )

def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
plt.show(fig1)

fig1.savefig('Pattern_Separation_BARS_MC_delete_'+str(Ndend)+'.eps', format = 'eps', dpi = 1200)
plt.close(fig1)


fig2, ax = plt.subplots()
ax.errorbar(output_MOSSY[0], output_MOSSY[1], yerr=output_MOSSY[2], color = 'green')
ax.errorbar(output_NO_MOSSY[0], output_NO_MOSSY[1], yerr=output_NO_MOSSY[2], color = 'lightgreen')


x_min = floor(min(min(output_MOSSY[0]), min(output_NO_MOSSY[0])))-5
x_max = ceil(max(max(output_MOSSY[0]), max(output_NO_MOSSY[0])))+5

y_min = floor(min(min(output_MOSSY[1]), min(output_NO_MOSSY[1])))-5
y_max = ceil(max(max(output_MOSSY[1]), max(output_NO_MOSSY[1])))+5

lims = [
    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
]


# now plot both limits against eachother
ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
#ax.set_aspect('equal')
plt.xlim(xmax = x_max, xmin = x_min)
plt.ylim(ymax = y_max, ymin = y_min)
plt.xticks([20])
plt.yticks(range(25, 81, 10))

plt.xticks(fontsize = FS-2)
plt.yticks(fontsize = FS-2)

xlabel('f1(input) (%)', fontsize = FS)
ylabel('f1(output) (%)', fontsize = FS)
legend( ('Control', 'Mossy loss') ,'lower right', fontsize = FS - 4)

plt.show(fig2)
fig2.savefig('Pattern_Separation_MC_delete_'+str(Ndend)+'.eps', format = 'eps', dpi = 1200)
plt.close(fig2)

fig3 = figure()

GC_MOSSY    = [x*100 for x in output_MOSSY[3]]
GC_NO_MOSSY = [x*100 for x in output_NO_MOSSY[3]]

plt.hist(GC_MOSSY, bins=20, histtype='stepfilled',
         weights=np.zeros_like(GC_MOSSY) + 1. / len(GC_MOSSY), color='green', alpha=0.8, label='Control')
plt.hist(GC_NO_MOSSY, bins=20, histtype='stepfilled',
         weights=np.zeros_like(GC_NO_MOSSY) + 1. / len(GC_NO_MOSSY), color='lightgreen', alpha=0.5, label='Mossy loss')

plt.xticks(fontsize = FS-2)
plt.yticks([x/100. for x in range(0,31,5)], fontsize = FS-2)

plt.xlabel("GC activity (%)", fontsize = FS)
plt.ylabel("Probability", fontsize = FS)
plt.legend( ('Control', 'Mossy loss') , fontsize = FS - 4)

plt.show(fig3)
fig3.savefig('GC_activity_MC_delete_'+str(Ndend)+'.eps', format = 'eps', bbox_inches='tight', dpi = 1200)
plt.close(fig3)

#==============================================================================
######################### S T A T I S T I C A L   T E S T S  ##################
#==============================================================================
# Statistical test for GC activity differences
sampleMOSSY = output_MOSSY[3]
sampleNO_MOSSY = output_NO_MOSSY[3]

meanMOSSY = mean(sampleMOSSY)
stdMOSSY  = std(sampleMOSSY)
seMOSSY   = stats.sem(sampleMOSSY)

meanNO_MOSSY = mean(sampleNO_MOSSY)
stdNO_MOSSY  = std(sampleNO_MOSSY)
seNO_MOSSY   = stats.sem(sampleNO_MOSSY)

###############################################################################
t_value = []
p_value = []

a = stats.ttest_ind(output_MOSSY[4][0], output_NO_MOSSY[4][0], axis=0, equal_var=False)
t_value.append(a[0])
p_value.append(a[1])

#Calculate the Wilcoxon signed-rank test.
#The Wilcoxon signed-rank test tests the null hypothesis that two related paired
#samples come from the same distribution. In particular, it tests whether the
#distribution of the differences x - y is symmetric about zero.
#It is a non-parametric version of the paired T-test.
T_value = []
p_value = []
a = stats.wilcoxon(output_MOSSY[4][0], output_NO_MOSSY[4][0], zero_method='wilcox', correction=True)
T_value.append(a[0])
p_value.append(a[1])

# Computes the Kolmogorov-Smirnov statistic on 2 samples.
# This is a two-sided test for the null hypothesis that 2 independent samples
# are drawn from the same continuous distribution.
a1 = stats.ks_2samp(sampleMOSSY, sampleNO_MOSSY)


mat1 = {'Control':sampleMOSSY, 'MC_delete': sampleNO_MOSSY}
np.save("GC_activity_MC_delete.npy", mat1)

mat2 = {'Control':output_MOSSY[4], 'MC_delete':output_NO_MOSSY[4]}
np.save("Similarity_out_MC_delete.npy", mat2)


rate1 = list()
for i in xrange(len(output_MOSSY[5])):
    rate1.append([2*x for x in output_MOSSY[5][i] if x != 0])
rate1 = list(itertools.chain(*rate1))
print mean(rate1)
print stats.sem(rate1)

rate2 = list()
for i in xrange(len(output_NO_MOSSY[5])):
    rate2.append([2*x for x in output_NO_MOSSY[5][i] if x != 0])
rate2 = list(itertools.chain(*rate2))
print mean(rate2)
print stats.sem(rate2)