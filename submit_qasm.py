import numpy as np
import random
import os
import glob

circuitPath = r"/home/haimeng/LocalProjects/IBM-PMME/Circuits"
device = "ibmq_bogota"
datastr_list = ['20220228', '20220301', '20220302']
QS = 'XplusState'
isDD = True
runtype = 'MeasMainqFreeLong' if not isDD else 'MeasMainqFreeLongXY4'

def shuffle_qasms(qasmlist):
    def get_time(filename):
        idgate = filename.split('/')[-1].split('_')[8].split('=')[-1]
        return idgate
    print(get_time(qasmlist[0]))
    times = list(set([get_time(x) for x in qasmlist])) #remove duplicates
    times.sort()
    print('# of time instances=%d'%(len(times)))
    random.shuffle(times) #randomize the timesteps
    print(times)
    sorted_qasm_list = []
    for time in times:
        sorted_qasm_list += [f for f in qasmlist if get_time(f) == time]
    return sorted_qasm_list

def get_qasms(datestr_list, runtype):
    runtypes_fullpath = []
    for datestr in datestr_list:
        allruns = os.scandir(circuitPath + '/' + device + '/' + datestr)
        runtypes_fullpath += [f.path for f in allruns if
                    f.is_dir() and QS in f.name and runtype+'_' in f.name]
    if len(set(runtypes_fullpath)) != 5:
        raise ValueError('total number of config incorrect')
    batchFiles = []
    batchFiles_flatten = []
    all_qasms = []
    # combine all qasms without measurement mitigation
    for runtype_fullpath in runtypes_fullpath:
        runtype = runtype_fullpath.split('/')[-1].split('_')[0]
        qasmlist = glob.glob(runtype_fullpath + '/%s*.qasm' % (runtype.split('_')[0]))
        if len(qasmlist) != 72:
            raise ValueError('incorrect circuit number in %s'%(runtype_fullpath))
        all_qasms += qasmlist
    if len(set(all_qasms)) != 5*72:
        raise ValueError('total number of qasms incorrect')
    for i in range(5):
        qasm_in_job = list(np.random.choice(all_qasms, size=72, replace=False))
        for q in qasm_in_job:
            all_qasms.remove(q)
        # add measurement calibration
        meas_qasm = glob.glob(runtype_fullpath + '/%s*.qasm' % ('MeasError'))
        qasm_in_job += meas_qasm
        batchFiles.append(qasm_in_job)
    for item in batchFiles:
        batchFiles_flatten += item
    if len(set(batchFiles_flatten)) != 5*72 + 2 or len(batchFiles_flatten) != 5*74:
        raise ValueError('total number of circuits incorrect')
    return batchFiles


def main():
    batchFiles = get_qasms(['20220228', '20220301','20220302'], runtype)
    return batchFiles


if __name__ == '__main__':
    main()