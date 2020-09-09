from datetime import date
import numpy as np
from qiskit.circuit import library
from qiskit.qasm import qasm
from math import pi
import os

rdigit = 5
prep_params = {"XplusState": (str(round(pi,rdigit)),str(0)),
               'XminusState':()}

meas_params = {"obsZ": None,
        "obsX": (str(round(pi,rdigit)),str(0)),
        "obsY": (str(0),str(round(3*pi/2,rdigit)))}

device_numQubit = {'ibmqx2':5,
                          'ibmq_armonk':1}

ddGateStr = {'X': 'u3(3.141592653589793,0,3.141592653589793)',
           'Y': 'u3(3.141592653589793,1.5707963267948966,1.5707963267948966)'}

def construct_exp_dict(runname,device,pstate,numIdGates,mbasis,circuitPath):
    """
    make the experiment configuration dict
    """
    today = date.today()
    datestr = today.strftime("%d%m%Y")
    exp_dict = {"runname": runname, "date": datestr, "device": device, "numIdGates": numIdGates,
                "pstate": pstate, "mbasis": mbasis, "circuitPath": circuitPath,
                "filename": runname + '_' + datestr + '_' + device + '_' + 'numIdGates=' + str(
                    numIdGates) + '_' + pstate + '_' + mbasis + '.qasm',
                "filepath": circuitPath + device + '/' + datestr + '/' + runname + '/'}
    return exp_dict


def construct_dd_exp_dict(runname,device,pstate,numDDrep,mbasis,circuitPath):
    """
    make the experiment configuration dict
    """
    today = date.today()
    datestr = today.strftime("%d%m%Y")
    exp_dict = {"runname": runname, "date": datestr, "device": device, "numDDrep": numDDrep,
                "pstate": pstate, "mbasis": mbasis, "circuitPath": circuitPath,
                "filename": runname + '_' + datestr + '_' + device + '_' + 'numXYXYseq=' + str(
                    numDDrep) + '_' + pstate + '_' + mbasis + '.qasm',
                "filepath": circuitPath + device + '/' + datestr + '/' + runname + '/'}
    return exp_dict


def write_qubit_index(numQubit):
    qstr = ''
    for i in range(numQubit):
        qstr += 'q[%s],'%(i)
    return qstr[:-1]

def write_state_tomo_qasms(**kwargs):
    if "exp_dict" in kwargs:
        exp_dict = kwargs.get("exp_dict")
    else:
        print("input dictionary name error")
    if not os.path.exists(exp_dict["filepath"]):
        os.makedirs(exp_dict["filepath"])
    numQubits = device_numQubit[exp_dict['device']]
    barrierStr = 'barrier %s;\n'%(write_qubit_index(numQubits))
    idStr = lambda q: 'id q[%s];\n'%(q)
    f = open(exp_dict["filepath"] + exp_dict["filename"], "w")
    # header for 5 qubit
    f.write("OPENQASM 2.0;\ninclude\"qelib1.inc\";\nqreg q[%s];\ncreg c[%s];\n"%(numQubits,numQubits))
    # state preparation
    prepGateStr = 'u2('+prep_params[exp_dict["pstate"]][0]+','+prep_params[exp_dict["pstate"]][1]+')'
    for i in range(numQubits):
        qubitStr = 'q[%d]'%(i)
        f.write(prepGateStr + ' ' + qubitStr +';\n')
    for d in range(exp_dict["numIdGates"]):
        f.write(barrierStr)
        for i in range(numQubits):
            f.write(idStr(i))
    # pre-measurement rotation
    if meas_params[exp_dict["mbasis"]] is not None:
        f.write(barrierStr)
        measGateStr = 'u2('+meas_params[exp_dict["mbasis"]][0]+','+meas_params[exp_dict["mbasis"]][1]+')'
        for i in range(numQubits):
            qubitStr = 'q[%d]' %(i)
            f.write(measGateStr + ' ' + qubitStr + ';\n')
    # measurement
    f.write(barrierStr)
    for i in range(numQubits):
        f.write('measure q[%d] -> c[%d];\n'%(i,i))
    f.close

def write_measurement_error_mitigation_qasm(**kwargs):
    if "exp_dict" in kwargs:
        exp_dict = kwargs.get("exp_dict")
    else:
        print("input dictionary name error")
    if not os.path.exists(exp_dict["filepath"]):
        os.makedirs(exp_dict["filepath"])
    numQubits = device_numQubit[exp_dict['device']]
    barrierStr = 'barrier %s;\n'%(write_qubit_index(numQubits))
    f = open(exp_dict["filepath"] + exp_dict["filename"], "w")
    # header for 5 qubit
    f.write("OPENQASM 2.0;\ninclude\"qelib1.inc\";\nqreg q[%s];\ncreg c[%s];\n"%(numQubits,numQubits))
    # z basis state prep

    # measurement basis post rotation
    # measure



def write_state_tomo_dd_qasms(**kwargs):
    if "exp_dict" in kwargs:
        exp_dict = kwargs.get("exp_dict")
    else:
        print("input dictionary name error")
    if not os.path.exists(exp_dict["filepath"]):
        os.makedirs(exp_dict["filepath"])
    numQubits = device_numQubit[exp_dict['device']]
    barrierStr = 'barrier %s;\n'%(write_qubit_index(numQubits))
    DDseq = ['X', 'Y', 'X', 'Y']
    f = open(exp_dict["filepath"] + exp_dict["filename"], "w")
    # header for 5 qubit
    f.write("OPENQASM 2.0;\ninclude\"qelib1.inc\";\nqreg q[%s];\ncreg c[%s];\n"%(numQubits,numQubits))
    # state preparation
    prepGateStr = 'u2('+prep_params[exp_dict["pstate"]][0]+','+prep_params[exp_dict["pstate"]][1]+')'
    for i in range(numQubits):
        qubitStr = 'q[%d]'%(i)
        f.write(prepGateStr + ' ' + qubitStr +';\n')
    for d in range(exp_dict["numDDrep"]):
        for i in range(len(DDseq)): #XYXY sequence, or XY$
            f.write(barrierStr)
            for j in range(numQubits):
                f.write(ddGateStr[DDseq[i]] + ' q[%s];\n'%(j))
    # pre-measurement rotation
    if meas_params[exp_dict["mbasis"]] is not None:
        f.write(barrierStr)
        measGateStr = 'u2('+meas_params[exp_dict["mbasis"]][0]+','+meas_params[exp_dict["mbasis"]][1]+')'
        for i in range(numQubits):
            qubitStr = 'q[%d]' %(i)
            f.write(measGateStr + ' ' + qubitStr + ';\n')
    # measurement
    f.write(barrierStr)
    for i in range(numQubits):
        f.write('measure q[%d] -> c[%d];\n'%(i,i))
    f.close

circuitPath = r"../Circuits/"
device = "ibmq_armonk"

runname = 'stateTomo_freeEvo'
measurement_basis = list(meas_params.keys())
# folder /device/date/state_tomography_freeEvo/
pstate = 'XplusState'

# # dense sampling rate
num_repetition = 48
num_complete = 0 # completed circuits
sampling_rate = 12
for r in range(num_complete,num_repetition):
    numIdGates = sampling_rate * r
    for mbasis in measurement_basis:
        dict0 = construct_exp_dict(runname=runname,device=device,pstate=pstate,mbasis=mbasis,circuitPath=circuitPath,numIdGates=numIdGates)
        write_state_tomo_qasms(exp_dict=dict0)

# # sparse sampling rate
num_repetition = 61
num_complete = 0 # completed circuits, starts with 0
sampling_rate = 24
offset = 600
for r in range(num_complete,num_repetition):
    numIdGates = sampling_rate * r + offset
    for mbasis in measurement_basis:
        dict0 = construct_exp_dict(runname=runname,device=device,pstate=pstate,mbasis=mbasis,circuitPath=circuitPath,numIdGates=numIdGates)
        write_state_tomo_qasms(exp_dict=dict0)

## DD circtuis
# runname = 'stateTomo_DD'
# measurement_basis = list(meas_params.keys())
# # folder /device/date/state_tomography_freeEvo/
# pstate = 'XplusState'
# num_repetition = 301
# num_complete = 0
# for r in range(num_complete,num_repetition):
#     for mbasis in measurement_basis:
#         dict0 = construct_dd_exp_dict(runname=runname,device=device,pstate=pstate,mbasis=mbasis,circuitPath=circuitPath,numDDrep=r)
#         write_state_tomo_dd_qasms(exp_dict=dict0)



# check a random qasm
# qasm_name = dict0["filepath"] + dict0["filename"]
# qasm_obj = qasm.Qasm(qasm_name)
# qasm_obj.parse()
#
# qasm_name2 = '/Users/Haimeng/Dropbox/IBM_PMME_data/Circuits/ibmq_armonk/11062020/stateTomo_freeEvo/' + 'stateTomo_freeEvo_11062020_ibmq_armonk_numIdGates=4_XplusState_obsX.qasm'
# qasm_obj2 = qasm.Qasm(qasm_name2)
# qasm_obj2.parse()