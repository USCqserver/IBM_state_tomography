"""
process tomogrpahy data in batch, test run in terminal
input: raw data directory
output:
save file in

"""
from tomography_fitting import *
TomoData = namedtuple('TomoData',['device','date','runname','qindex','pstate','spstate','runNo','jobNo'])
batch_list = [TomoData('ibmq_athens','20210302','MeasMainqFreeLong',0,'XplusState','ZplusState','run1',p) for p in range(1,5)]
batch_list += [TomoData('ibmq_athens','20210302','MeasMainqFreeShort',0,'XplusState','ZplusState','run1',p) for p in range(5)]
# batch_list = [TomoData('ibmq_athens','20210302','MeasMainqFreeLong',0,'XplusState','ZminusState','run1',p) for p in range(5)]
# batch_list += [TomoData('ibmq_athens','20210302','MeasMainqFreeShort',0,'XplusState','ZminusState','run1',p) for p in range(5)]
# check if data folder exists
for b in batch_list:
    exp0 = ExpObject(runname='%s_Q%d_%s_QS_%s_%d' % (b.runname,b.qindex,b.pstate, b.spstate,b.jobNo), datestr=b.date, device=b.device,
                     pstate=b.pstate,
                     runNo=b.runNo)
    datapath = r"/home/haimeng/LocalProjects/IBM-PMME/Data/raw/" + exp0.device + "/" + exp0.datestr + "/" + exp0.runname  + "/" + exp0.runNo + "/"
    if not glob.glob(datapath):
        raise ValueError('%s does not exists'%(datapath))

for b in batch_list:
    exp0 = ExpObject(runname='%s_Q%d_%s_QS_%s_%d'%(b.runname,b.qindex,b.pstate,b.spstate,b.jobNo), datestr=b.date, device=b.device, pstate=b.pstate,
                     runNo=b.runNo)
    datapath = r"/home/haimeng/LocalProjects/IBM-PMME/Data/raw/" + exp0.device + "/" + exp0.datestr + "/" + exp0.runname + "/" + exp0.runNo +"/"
    runname = glob.glob(datapath)[0].split('/')[-3]
    runtype = runname.split('_')[0]
    # store data to
    filepath = r"/home/haimeng/LocalProjects/IBM-PMME/Analysis/" + exp0.device + "/" + exp0.datestr + '/' +runtype + '/'
    tomo_fit = StateTomographyFit()
    save = True
    Nboots = 100
    tomo_fit.fit_state_from_run(datapath,exp0,qubiti=[b.qindex],polar=False,strid=8,readout_error_mitigation=readout_error_mitigate)
    tomo_fit.blochInPolarCoordinate(save=save,show=False,filepath=filepath,readout_error_mitigate=readout_error_mitigate)
    # change where x-axis stops using keyword argument 'idMax';
    # change x-axis resolution using argument 'spacing', sampling frequency equals 1 sample per spacing*4 Id gates
    tomo_fit.assign_errorbars_from_run(exp0,qubiti=[b.qindex],Nboots=Nboots,readout_error_mitigation=readout_error_mitigate,bootstrap_method=bootstrap_method)
    tomo_fit.plot_bloch_vector(exp0,save=save, show=False,filepath=filepath, spacing=1, polar=False,
                               readout_error_mitigation=readout_error_mitigate)
    csvfile = tomo_fit.saveBloch2csv(filepath, bootstrap=True,polar=True, readout_error_mitigation=readout_error_mitigate)