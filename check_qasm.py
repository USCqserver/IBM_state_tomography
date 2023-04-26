"""
input exp dir of 74 circtuis under one config
"""
import glob
from qiskit import QuantumCircuit
circuitPath = '/home/haimeng/LocalProjects/IBM-PMME/Circuits'
device = 'ibmq_bogota'
datestr1 = '20220301'
datestr2 = '20220228'
datestr3 = '20220302'
runs = glob.glob('/'.join([circuitPath, device, datestr1, '*'])) + \
    glob.glob('/'.join([circuitPath, device, datestr2,'*'])) + \
    glob.glob('/'.join([circuitPath, device, datestr3,'*']))
runtype = 'MeasMainqFreeLong_Q0_randomState_QS4_XplusState'
batchFiles = []
qasmpath =runs[1]
qasmlist = glob.glob(qasmpath + '/%s*.qasm'%(runtype.split('_')[0]))
# sort by id gates
qasmlist.sort(key=lambda x: int(x.split('/')[-1].split('_')[8].split('=')[-1]))
print(len(qasmlist))
qname =qasmlist[4]
print(qname.split('/')[-1])
f = QuantumCircuit().from_qasm_file(qname)
f.draw()
# prep_key =
# meas_key =
# count each operation
# number of id
# number of dd

# count circuit depth
