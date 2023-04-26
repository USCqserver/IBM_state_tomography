QS = 'ZplusState'
import glob
import os, shutil
from pathlib import Path
runname = 'MeasMainqFreeLong'
runNo = 'run1'
tetra_states = ['tetra1', 'tetra2', 'tetra3', 'randomState']
datadir = '/home/haimeng/LocalProjects/IBM-PMME/Data/raw/'
device = 'ibmq_bogota'
datestr_list = ['20220228', '20220302']
for p in tetra_states[-1:]:
    directories = []
    for i in range(5):
        directories += glob.glob(datadir + device + '/' + '*' + '/' + '%s_Q0_%s_QS4_%s_%d'%(runname, p, QS,i) + '/'+runNo + '/')
    if len(directories)!=5:
        raise ValueError('not enough job folder found for state %s'%(p))
    dest = datadir + device + '/' + '20220307' + '/' + '%s_Q0_%s_QS4_%s_all'%(runname, p, QS) + '/' + runNo +'/'
    if not os. path.exists(dest):
        Path(dest).mkdir(parents=True, exist_ok=True)
    all_files = []
    for dir in directories:
        all_files += glob.glob(dir + '*.txt')
    # fileNum = len([name for name in os.listdir(dest) if os.path.isfile(os.path.join(dest,name))])
    if len(all_files) != 74:
        raise ValueError('not enough files')
    else:
        for file in all_files:
            shutil.copy(file, dest)