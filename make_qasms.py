from datetime import date
import numpy as np
from qiskit.circuit import library
from qiskit.qasm import qasm
from math import pi
import os

rdigit = 5
prep_params = {"XplusState": ('u2',str(0),str(round(pi,rdigit))),
               'XminusState':('u2',str(round(pi,rdigit)),str(0)),
               'YplusState':('u2',str(round(pi/2,rdigit)),str(round(pi/2,rdigit))),
               'YminusState':('u2',str(round(-pi/2,rdigit)),str(round(-pi/2,rdigit))),
               'ZplusState':('id',''),
               'ZminusState':('u3',str(round(pi,rdigit)),str(-round(pi/2,rdigit)),str(round(pi/2,rdigit)))}
# measurement gates are all u2 gate
meas_params = {"obsZ": None,
        "obsX": (str(0),str(round(pi,rdigit))),
        "obsY": (str(0),str(round(pi/2,rdigit)))}

device_numQubit = {'ibmqx2':5,
                          'ibmq_armonk':1}

ddGateStr = {'X': 'u3(3.141592653589793,0,3.141592653589793)',
           'Y': 'u3(3.141592653589793,1.5707963267948966,1.5707963267948966)'}

def construct_exp_dict(device,pstate,mbasis,circuitPath,**kwargs):
    """
    make the experiment configuration dict for free evolution
    """
    today = date.today()
    datestr = today.strftime("%d%m%Y")
    if 'runname' in kwargs:
        runname = kwargs.get('runname')
    elif 'mainq' in kwargs and isinstance(pstate,dict):
        runname = '_'.join(['Free','Q%d'%(int(kwargs.get('mainq'))),pstate['main'],'QS',pstate['spec']])
    elif isinstance(pstate,str):
        runname = 'Free' + '_' + pstate
    exp_dict = {"runname": runname, "date": datestr, "device": device,
                 "circuitPath": circuitPath,'pstate': pstate, 'mbasis': mbasis,
                "filepath": circuitPath + device + '/' + datestr + '/' + runname + '/'}
    if 'mainq' in kwargs:
        exp_dict['mainq'] = int(kwargs.get('mainq'))
    if 'numIdGates' in kwargs:
        exp_dict['numIdGates'] = kwargs.get('numIdGates')
        exp_dict["filename"]= runname + '_' + datestr + '_' + device + '_' + 'numIdGates=' + str(
            exp_dict['numIdGates']) + '_' + exp_dict['mbasis'] + '.qasm'
    else:
        exp_dict["filename"] = runname + '_' + datestr + '_' + device + '_' + exp_dict['pstate'] + '_' + exp_dict['mbasis'] + '.qasm'
    return exp_dict


def construct_dd_exp_dict(device,pstate,numDDrep,mbasis,circuitPath,**kwargs):
    """
    make the experiment configuration dict
    """
    today = date.today()
    datestr = today.strftime("%d%m%Y")
    if 'runname' in kwargs:
        runname = kwargs.get('runname')
    elif 'mainq' in kwargs and isinstance(pstate,dict):
        runname = '_'.join(['DD','Q%d'%(int(kwargs.get('mainq'))),pstate['main'],'QS',pstate['spec']])
    elif isinstance(pstate,str):
        runname = 'DD' + '_' + pstate
    exp_dict = {"runname": runname, "date": datestr, "device": device, "numDDrep": numDDrep,
                "pstate": pstate, "mbasis": mbasis, "circuitPath": circuitPath,
                "filename": runname + '_' + datestr + '_' + device + '_' + 'numXYXYseq=' + str(
                    numDDrep) + '_' + mbasis + '.qasm',
                "filepath": circuitPath + device + '/' + datestr + '/' + runname + '/'}
    if 'mainq' in kwargs:
        exp_dict['mainq'] = kwargs.get('mainq')
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
    if prep_params[exp_dict["pstate"]] is not None:
        gate = prep_params[exp_dict["pstate"]][0]
        if gate == 'u2':
            prepGateStr = 'u2('+prep_params[exp_dict["pstate"]][1]+','+prep_params[exp_dict["pstate"]][2]+')'
        elif gate =='u3':
            prepGateStr = 'u3(' + prep_params[exp_dict["pstate"]][1] + ',' + prep_params[exp_dict["pstate"]][2]+ ',' + prep_params[exp_dict["pstate"]][3] + ')'
        elif gate=='id':
            prepGateStr = 'id'
        else:
            raise NameError('gate name does not match names in gate library')
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

def write_measurement_error_mitigation_qasm(cal_state,cal_basis,**kwargs):
    if "exp_dict" in kwargs:
        exp_dict = kwargs.get("exp_dict")
    else:
        print("input dictionary name error")
    if not os.path.exists(exp_dict["filepath"]):
        os.makedirs(exp_dict["filepath"])
    numQubits = device_numQubit[exp_dict['device']]
    barrierStr = 'barrier %s;\n'%(write_qubit_index(numQubits))
    filename = exp_dict['filename'].replace(exp_dict['runname'],'MeasError_Mitigate')
    f = open(exp_dict["filepath"] + filename, "w")
    # header for 5 qubit
    f.write("OPENQASM 2.0;\ninclude\"qelib1.inc\";\nqreg q[%s];\ncreg c[%s];\n"%(numQubits,numQubits))
    # z basis state prep
    # state preparation
    if prep_params[cal_state] is not None:
        gate = prep_params[cal_state][0]
        if gate == 'u2':
            prepGateStr = 'u2(' + prep_params[cal_state][1] + ',' + prep_params[cal_state][2] + ')'
        elif gate == 'u3':
            prepGateStr = 'u3(' + prep_params[cal_state][1] + ',' + prep_params[cal_state][2] + ',' + \
                          prep_params[exp_dict["pstate"]][3] + ')'
        elif gate == 'id':
            prepGateStr = 'id'
        i = exp_dict['mainq']
        qubitStr = 'q[%d]' % (i)
        f.write(prepGateStr + ' ' + qubitStr + ';\n')
    # measurement basis post rotation
    if meas_params[cal_basis] is not None:
        f.write(barrierStr)
        measGateStr = 'u2('+meas_params[cal_basis][0]+','+meas_params[cal_basis][1]+')'
        qubitStr = 'q[%d]' %(i)
        f.write(measGateStr + ' ' + qubitStr + ';\n')
    # measurement
    f.write(barrierStr)
    # for i in range(numQubits):
    f.write('measure q[%d] -> c[%d];\n'%(i,i))
    f.close


def write_state_tomo_dd_qasms(**kwargs):
    """
    DD on all qubits
    """
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
    gate = prep_params[exp_dict["pstate"]][0]
    if gate == 'u2':
        prepGateStr = 'u2('+prep_params[exp_dict["pstate"]][1]+','+prep_params[exp_dict["pstate"]][2]+')'
    elif gate == 'u3':
        prepGateStr = 'u3(' + prep_params[exp_dict["pstate"]][1] + ',' + prep_params[exp_dict["pstate"]][2] + ','+ prep_params[exp_dict["pstate"]][3] + ')'
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
def write_state_prep(f,exp_dict):
    """
    pstate is a string or a dictionary of strings
    """
    numQubits = numQubits = device_numQubit[exp_dict['device']]
    if 'mainq' in exp_dict:
        mainq = exp_dict['mainq']
    for i in range(numQubits):
        qubitStr = 'q[%d]' % (i)
        # get the gate type
        if i == mainq:
            gatename = exp_dict["pstate"]['main']
            gate = prep_params[gatename][0]
        else:
            gatename = exp_dict["pstate"]['spec']
        gate = prep_params[gatename][0]
        if gate == 'u2':
            prepGateStr = 'u2(' + prep_params[gatename][1] + ',' + prep_params[gatename][
                2] + ')'
        elif gate == 'u3':
            prepGateStr = 'u3(' + prep_params[gatename][1] + ',' + prep_params[gatename][
                2] + ',' + prep_params[gatename][3] + ')'
        f.write(prepGateStr + ' ' + qubitStr + ';\n')

def write_state_tomo_free_qasms(**kwargs):
    """
    advance version of free evolution with option to prepare spectator qubits with user-specified initial states
    exp_dict[pstate] is a dictionary contains initial states for main qubit and spectator qubits
    """
    if "exp_dict" in kwargs:
        exp_dict = kwargs.get("exp_dict")
    else:
        print("input dictionary name error")
    if not os.path.exists(exp_dict["filepath"]):
        os.makedirs(exp_dict["filepath"])
    mainq = int(exp_dict['mainq']) # get the mainq index
    numQubits = device_numQubit[exp_dict['device']]
    barrierStr = 'barrier %s;\n'%(write_qubit_index(numQubits))
    idStr = lambda q: 'id q[%s];\n'%(q)
    f = open(exp_dict["filepath"] + exp_dict["filename"], "w")
    # header for 5 qubit
    f.write("OPENQASM 2.0;\ninclude\"qelib1.inc\";\nqreg q[%s];\ncreg c[%s];\n"%(numQubits,numQubits))
    # state preparation
    write_state_prep(f,exp_dict)
    # free evolution
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
            if i == mainq:
                f.write(measGateStr + ' ' + qubitStr + ';\n')
    # measurement
    f.write(barrierStr)
    # for i in range(numQubits):
    i = exp_dict['mainq']
    f.write('measure q[%d] -> c[%d];\n'%(i,i))
    f.close

def write_state_tomo_dd_on_spectator_qasms(**kwargs):
    """
    DD on all qubits
    """
    if "exp_dict" in kwargs:
        exp_dict = kwargs.get("exp_dict")
    else:
        print("input dictionary name error")
    if not os.path.exists(exp_dict["filepath"]):
        os.makedirs(exp_dict["filepath"])
    numQubits = device_numQubit[exp_dict['device']]
    mainq = int(exp_dict['mainq'])  # get the mainq index
    barrierStr = 'barrier %s;\n'%(write_qubit_index(numQubits))
    DDseq = ['X', 'Y', 'X', 'Y']
    f = open(exp_dict["filepath"] + exp_dict["filename"], "w")
    # header for 5 qubit
    f.write("OPENQASM 2.0;\ninclude\"qelib1.inc\";\nqreg q[%s];\ncreg c[%s];\n"%(numQubits,numQubits))
    # state preparation
    write_state_prep(f,exp_dict)
    for d in range(exp_dict["numDDrep"]):
        for i in range(len(DDseq)): #XYXY sequence, or XY$
            f.write(barrierStr)
            for j in range(numQubits):
                if j != mainq:
                    f.write(ddGateStr[DDseq[i]] + ' q[%s];\n'%(j))
    # pre-measurement rotation
    if meas_params[exp_dict["mbasis"]] is not None:
        f.write(barrierStr)
        measGateStr = 'u2('+meas_params[exp_dict["mbasis"]][0]+','+meas_params[exp_dict["mbasis"]][1]+')'
        for i in range(numQubits):
            qubitStr = 'q[%d]' %(i)
            if i == mainq:
                f.write(measGateStr + ' ' + qubitStr + ';\n')
    # measurement
    f.write(barrierStr)
    # for i in range(numQubits):
    i = exp_dict['mainq']
    f.write('measure q[%d] -> c[%d];\n'%(i,i))
    f.close

circuitPath = r"../Circuits/"
device = "ibmqx2"

runname = 'stateTomo_freeEvo'
measurement_basis = list(meas_params.keys())
# folder /device/date/state_tomography_freeEvo/
# pstates = ['XminusState', 'YminusState','YplusState']
pstates = ['XplusState']
def write_free_and_dd():
    for pstate in pstates:
        runname = 'stateTomo_freeEvo'
        runname = runname + '_' + pstate
        # measurement mitigation circuits
        for p in ['ZplusState','ZminusState']:
            dict0 = construct_exp_dict(runname=runname,device=device,pstate=p,mbasis='obsZ',circuitPath=circuitPath)
            write_measurement_error_mitigation_qasm(p,'obsZ',exp_dict=dict0)
        # # dense sampling rate
        num_repetition = 24
        num_complete = 0 # completed circuits
        sampling_rate = 24
        for r in range(num_complete,num_repetition):
            numIdGates = sampling_rate * r
            for mbasis in measurement_basis:
                dict0 = construct_exp_dict(runname=runname,device=device,pstate=pstate,mbasis=mbasis,circuitPath=circuitPath,numIdGates=numIdGates)
                write_state_tomo_qasms(exp_dict=dict0)

    # # # sparse sampling rate
    # num_repetition = 61
    # num_complete = 0 # completed circuits, starts with 0
    # sampling_rate = 24
    # offset = 600
    # for r in range(num_complete,num_repetition):
    #     numIdGates = sampling_rate * r + offset
    #     for mbasis in measurement_basis:
    #         dict0 = construct_exp_dict(runname=runname,device=device,pstate=pstate,mbasis=mbasis,circuitPath=circuitPath,numIdGates=numIdGates)
    #         write_state_tomo_qasms(exp_dict=dict0)

        # DD circtuis
        runname = 'stateTomo_DD' + '_' + pstate
        for p in ['ZplusState','ZminusState']:
            dict0 = construct_exp_dict(runname=runname,device=device,pstate=p,mbasis='obsZ',circuitPath=circuitPath)
            write_measurement_error_mitigation_qasm(p,'obsZ',exp_dict=dict0)
        for r in range(num_complete,num_repetition):
            numDDgates = sampling_rate * r
            for mbasis in measurement_basis:
                dict0 = construct_dd_exp_dict(runname=runname,device=device,pstate=pstate,mbasis=mbasis,circuitPath=circuitPath,numDDrep=numDDgates)
                write_state_tomo_dd_qasms(exp_dict=dict0)

def write_free_and_dd_w_spectators(qindex,pstate):
        # # dense sampling rate
        num_repetition = 24
        num_complete = 0  # completed circuits
        sampling_rate = 24
        for r in range(num_complete, num_repetition):
            numIdGates = sampling_rate * r
            for mbasis in measurement_basis:
                dict0 = construct_exp_dict(mainq=qindex, device=device, pstate=pstate, mbasis=mbasis,
                                           circuitPath=circuitPath, numIdGates=numIdGates)
                write_state_tomo_free_qasms(exp_dict=dict0)
        # measurement mitigation circuits
        for p in ['ZplusState', 'ZminusState']:
            dict0 = construct_exp_dict(mainq=qindex,runname=dict0['runname'], device=device, pstate=p, mbasis='obsZ',
                                       circuitPath=circuitPath)
            write_measurement_error_mitigation_qasm(p, 'obsZ', exp_dict=dict0)
        # DD circtuis
        for r in range(num_complete, num_repetition):
            numDDgates = sampling_rate * r
            for mbasis in measurement_basis:
                dict0 = construct_dd_exp_dict(mainq=qindex,device=device, pstate=pstate, mbasis=mbasis,
                                              circuitPath=circuitPath, numDDrep=numDDgates)
                write_state_tomo_dd_on_spectator_qasms(exp_dict=dict0)
        for p in ['ZplusState', 'ZminusState']:
            dict0 = construct_exp_dict(mainq=qindex,runname=dict0['runname'], device=device, pstate=p, mbasis='obsZ',
                                       circuitPath=circuitPath)
            write_measurement_error_mitigation_qasm(p, 'obsZ', exp_dict=dict0)
for p in prep_params.keys():
    pstates = {'main': 'XplusState',
               'spec': p}
    write_free_and_dd_w_spectators(qindex=0,pstate=pstates)
# check a random qasm
# qasm_name = dict0["filepath"] + dict0["filename"]
# qasm_obj = qasm.Qasm(qasm_name)
# qasm_obj.parse()
# #
# qasm_name2 = '../Circuits/ibmqx2/09112020/stateTomo_freeEvo_XplusState/' + 'stateTomo_freeEvo_XplusState_09112020_ibmqx2_numIdGates=48_XplusState_obsX.qasm'
# qasm_obj2 = qasm.Qasm(qasm_name2)
# qasm_obj2.parse()