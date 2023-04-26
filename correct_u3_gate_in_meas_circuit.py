import numpy as np
import random
import os
import glob

circuitPath = r"/home/haimeng/LocalProjects/IBM-PMME/Circuits"
device = "ibmq_bogota"
datestr_list = ['20220228', '20220301', '20220302']
for datestr in datestr_list:
    allruns = os.scandir(circuitPath + '/' + device + '/' + datestr)
    for p in allruns:
        q = glob.glob(p.path + '/MeasError_Mitigate_*_ibmq_bogota_ZminusState_obsZ.qasm')
        f = open(q[0],'r+')
        qasm = f.read()
        f.truncate(0)
        f.seek(0)
        commands = qasm.split(';')
        if 'u3' in qasm:
            for i in range(len(commands)):
                if 'u3' in commands[i]:
                    commands[i] = 'rz(-pi) q[0];\nsx q[0];\nsx q[0];\nrz(pi) q[0]'
            qasm = ';'.join(commands) + ';'
            f.write(qasm)
        else:
            qasm = ';'.join(commands[:-1])+';'
            f.write(qasm)
        f.close()

