from measurement_error_mitigation import *
from tomography_fitting import *

def generate_mem_data_from_raw(data,datadir):
    if 'Spec' in data.runname:
        runname = '%s_Q%d_%s_QS%d_%s_%d' % (data.runname, data.qindex, data.pstate, 4-data.jobNo,data.spstate, data.jobNo)
    else:
        runname = '%s_Q%d_%s_QS_%s_%d' % (data.runname, data.qindex, data.pstate, data.spstate, data.jobNo)
    exp0 = ExpObject(runname=runname, datestr=data.date,
                     device=data.device,
                     pstate=data.pstate,
                     runNo=data.runNo)
    datapath = datadir + exp0.device + "/" + exp0.datestr + "/" + exp0.runname + "/" + exp0.runNo + "/"
    if not glob.glob(datapath):
        raise ValueError('%s does not exists'%(datapath))
    mem_files = glob.glob(datapath+'MeasError_Mitigate*ZplusState*.txt') + glob.glob(datapath+'MeasError_Mitigate*ZminusState*.txt')
    mem_data = [countWithBinaryKey(f) for f in mem_files]
    mem_matrix = mem_matrix_from_data(mem_data)
    tomo_files = glob.glob(datapath+'%s*.txt'%data.runname)
    for file in tomo_files:
        data = countWithBinaryKey(file)
        new_dict = mem_on_data(data,mem_matrix)
        # save dict
        save_path = file.replace('raw','mem')
        save_dir = datapath.replace('raw','mem')
        if not os.path.exists(save_dir):
            Path(save_dir).mkdir(parents=True)
        lines = ['dummy line','dummy line'] + ['"%s"'%(key)+',%d'%(round(item)) for key,item in new_dict.items()]
        with open(save_path,'w') as f:
            for line in lines:
                f.write(line + '\n')

def main():
    for i in range(1,4):
        runNo = 'run%d'%(i)
        batch_list0 = [TomoData('ibmq_athens', '20210625', 'MeasMainqFreeLong', 0, p, s, runNo, i) for s in spec_pstates[1:]
                   for i,p in enumerate(tetra_pstates)]
        # batch_list0 = [TomoData('ibmq_athens', '20210625', 'NonMarkMeasMainqFreeLong', 0, p, s, runNo, i) for s in spec_pstates
        #            for i,p in enumerate(nonmark_pstates)]
        # batch_list0 = [TomoData('ibmq_athens', '20210625', 'SpecMeasMainqFreeLong', 0, 'tetra0', s, runNo, j) for s in spec_pstates
        #            for j in range(5)]
        datadir = r"/home/haimeng/LocalProjects/IBM-PMME/Data/raw/"
        for b in batch_list0:
            generate_mem_data_from_raw(b,datadir)

if __name__ == "__main__":
    main()
