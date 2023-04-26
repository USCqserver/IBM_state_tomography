datestr_list = ['20220307', '20220308']
from measurement_error_mitigation import *
from tomography_fitting import *

QS = 'ZplusState'
device = 'ibmq_bogota'
datestr = '20220307'
datadir = r"/home/haimeng/LocalProjects/IBM-PMME/Data/raw/"
runname = 'MeasMainqFreeLong'
def generate_mem_data_from_raw(index):
    datapath = datadir + device + "/" + datestr + "/" + '%s_QS=%s_job%d'%(runname, QS, index) + "/"
    if not glob.glob(datapath):
        raise ValueError('%s does not exists'%(datapath))
    mem_files = glob.glob(datapath+'MeasError_Mitigate*ZplusState*.txt') + glob.glob(datapath+'MeasError_Mitigate*ZminusState*.txt')
    mem_data = [countWithBinaryKey(f) for f in mem_files]
    mem_matrix = mem_matrix_from_data(mem_data)
    tomo_files = glob.glob(datapath+'%s*.txt'%runname)
    for file in tomo_files:
        data = countWithBinaryKey(file)
        new_dict = mem_on_data(data,mem_matrix)
        # save dict
        save_path = file.replace('raw','mem') # this one contains the filename
        save_path = save_path.split('/')
        save_dir = datapath.replace('raw','mem')
        if not os.path.exists(save_dir):
            Path(save_dir).mkdir(parents=True)
        lines = ['dummy line','dummy line'] + ['"%s"'%(key)+',%d'%(round(item)) for key,item in new_dict.items()]
        with open(save_path,'w') as f:
            for line in lines:
                f.write(line + '\n')


for datestr in datestr_list:
    # there are 5 jobs in one batch
    batch_list = [TomoData('ibmq_bogota', datestr, 'MeasMainqFreeLong', 0, p, QS, 'run1', i)
                   for i, p in enumerate(tetra_pstates)]
