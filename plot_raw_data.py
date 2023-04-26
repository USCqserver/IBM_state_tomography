# plot raw data
"""
input:
experiment configuration and data folder
output
id number and px, py, pz
"""
import glob
from tomography_fitting import make_tomography_hist, ExpObject, TomoData
from tomography_fitting import spec_pstates, tetra_pstates, nonmark_pstates
import numpy as np
import matplotlib.pyplot as plt

data_root = '../Data/raw/ibmq_athens/20210701/'
# for each file read the histgram
strid = 8
runNo = 'run2'
Spec_pstates = ['tetra0'] * 5
def plot_probabilities(ax,tlist,prob):
    mbasis = ['obsZ', 'obsX', 'obsY']
    for i in range(prob.shape[1]):
        ax.plot(tlist,prob[:,i],'-x',label = mbasis[i])
    return ax
fig1, axes = plt.subplots(3,5)
for i,QS in enumerate(spec_pstates):
    for jobid,Q0 in enumerate(Spec_pstates):
        runname = 'SpecRevMeasMainqFreeLong'
        # runname = 'MeasMainqFreeLong'
        txtfiles = glob.glob(data_root + '%s_Q0_%s_QS%d_%s_*/%s/%s*.txt' % (runname, Q0,4-jobid, QS, runNo,runname))
        # txtfiles = glob.glob(data_root + '%s_Q0_%s_QS_%s_*/%s/%s*.txt' % (runname, Q0,QS, runNo,runname))
        job = TomoData('ibmq_athens', '20210701', runname, 0, Q0, QS, runNo, jobid)
        exp0 = ExpObject(runname='%s_Q%d_%s_QS%d_%s' % (job.runname, job.qindex, job.pstate,4-i,job.spstate),
                         datestr=job.date, device=job.device, pstate=job.pstate,
                         

                         runNo=job.runNo)
        idlist = sorted(list(set([int(f.split('/')[-1].split('_')[strid].split('=')[1]) for f in txtfiles])))
        histograms = np.array(
            [make_tomography_hist(exp0, numIdGates=id, qubiti=[0], shots=8192) for id in idlist]) / 8192

        axes[i,jobid] = plot_probabilities(axes[i,jobid], idlist, histograms[:,:,0])
axes[0,0].legend()
fig1.show()
savedir = r'../Analysis/ibmq_athens/20210701/plots/'
fig1.savefig(savedir + runname + '_' + runNo)
