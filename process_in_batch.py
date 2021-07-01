"""
process tomogrpahy data in batch, test run in terminal
input: raw data directory
output:
save file in
environment:
"""
from tomography_fitting import *
import time
startTime = time.time()
save = True
Nboots = 100
bootsMehtod = 'bayesian'


# tetra_pstates = ['tetra1']
spec_pstates = [ 'XplusState']
tetra_pstates = ['tetra3', 'tetra2', 'tetra1', 'tetra0', 'randomState']
nonmark_pstates = ['XplusState', 'XminusState', 'tetra0', 'YplusState', 'YminusState']
batch_list0 = [TomoData('ibmq_athens','20210625','MeasMainqFreeLong',0,p,'XplusState','run1',i) for i,p in enumerate(tetra_pstates)]
batch_list1 = [TomoData('ibmq_athens','20210625','NonMarkMeasMainqFreeLong',0, p ,s ,'run1',i) for s in spec_pstates for i,p in enumerate(nonmark_pstates)]
batch_list2 = [TomoData('ibmq_athens','20210625','SpecMeasMainqFreeLong',0,'tetra0', s ,'run1',i) for s in spec_pstates for i,p in enumerate(tetra_pstates)]
# batch_list = [TomoData('ibmq_athens','20210302','MeasMainqFreeLong',0,'XplusState','ZplusState','run1',p) for p in range(5)]
# batch_list += [TomoData('ibmq_athens','20210302','MeasMainqFreeShort',0,'XplusState','ZplusState','run1',p) for p in range(5)]
# batch_list = [TomoData('ibmq_athens','20210302','MeasMainqFreeLong',0,'XplusState','ZminusState','run1',p) for p in range(5)]
# batch_list += [TomoData('ibmq_athens','20210302','MeasMainqFreeShort',0,'XplusState','ZminusState','run1',p) for p in range(5)]
# check if data folder exists
batch_list =  batch_list1[2:]
print('number of batches = %d'%(len(batch_list)))
for b in batch_list:
    if 'Spec' in b.runname:
        exp0 = ExpObject(runname='%s_Q%d_%s_QS%d_%s_%d'%(b.runname,b.qindex,b.pstate,4-int(b.jobNo),b.spstate, b.jobNo), datestr=b.date, device=b.device, pstate=b.pstate,
                     runNo=b.runNo)
    else:
        exp0 = ExpObject(runname='%s_Q%d_%s_QS_%s_%d' % (b.runname,b.qindex,b.pstate, b.spstate,b.jobNo), datestr=b.date, device=b.device,
                     pstate=b.pstate,
                     runNo=b.runNo)
    datapath = r"/home/haimeng/LocalProjects/IBM-PMME/Data/raw/" + exp0.device + "/" + exp0.datestr + "/" + exp0.runname + "/" + exp0.runNo + "/"
    if not glob.glob(datapath):
        raise ValueError('%s does not exists'%(datapath))

for b in batch_list:
    if 'Spec' in b.runname:
        exp0 = ExpObject(runname='%s_Q%d_%s_QS%d_%s_%d'%(b.runname,b.qindex,b.pstate,4-int(b.jobNo),b.spstate, b.jobNo), datestr=b.date, device=b.device, pstate=b.pstate,
                     runNo=b.runNo)
    else:
        exp0 = ExpObject(runname='%s_Q%d_%s_QS_%s_%d'%(b.runname,b.qindex,b.pstate,b.spstate,b.jobNo), datestr=b.date, device=b.device, pstate=b.pstate,
                     runNo=b.runNo)
    datapath = r"/home/haimeng/LocalProjects/IBM-PMME/Data/raw/" + exp0.device + "/" + exp0.datestr + "/" + exp0.runname + "/" + exp0.runNo +"/"
    runname = glob.glob(datapath)[0].split('/')[-3]
    fittype = 'method_%s_mem_%s'%(method,readout_error_mitigate)
    # store data to
    filepath = r"/home/haimeng/LocalProjects/IBM-PMME/Analysis/" + exp0.device + "/" + exp0.datestr + '/' +runname + '/' + fittype + '/'
    tomo_fit = StateTomographyFit(exp0,datapath,qubiti=[b.qindex],strid=8)
    tomo_fit.fit_state_from_run(polar=False,readout_error_mitigation=readout_error_mitigate)
    tomo_fit.blochInPolarCoordinate(save=save,show=False,filepath=filepath,readout_error_mitigate=readout_error_mitigate)
    # change where x-axis stops using keyword argument 'idMax';
    # change x-axis resolution using argument 'spacing', sampling frequency equals 1 sample per spacing*4 Id gates
    generate_bootstrapped_tomography_dataset(tomo_fit,filepath=filepath+'%s_bootstrapped/'%(bootsMehtod),Nboots=Nboots,readout_error_mitigation=readout_error_mitigate,bootstrap_method=bootsMehtod,save=save,filename=runname)
    tomo_fit.plot_bloch_vector(exp0,save=save, show=False,filepath=filepath, spacing=1, polar=False,
                               readout_error_mitigation=readout_error_mitigate)
    csvfile = tomo_fit.saveBloch2csv(filepath, readout_error_mitigation=readout_error_mitigate)

excutionTime = time.time()-startTime
print('%s batches of data, %d bootstrap samples, %s bootstrap method takes: %.2f min'%(len(batch_list), Nboots, bootsMehtod, excutionTime/60))
