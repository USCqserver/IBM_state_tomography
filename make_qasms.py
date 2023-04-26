"""
starting Feb 28, 2022 run in qiskit 0.23.6
before running in 0.19.1
"""
from datetime import date
import numpy as np
from qiskit.circuit import library
from qiskit import transpile
from qiskit.circuit.quantumcircuit import QuantumCircuit
from qiskit.qasm import qasm
from qiskit import IBMQ
from qiskit.quantum_info import operators, Operator
from qiskit.extensions import unitary
from math import pi
import os
from utility import cartesian_to_polar, myStr2float, paulix, pauliy, pauliz
import re
rdigit = 5
# measurement gates are all u2 gate
meas_params = {"obsZ": None,
        "obsX": (str(0),str(round(pi,rdigit))),
        "obsY": (str(0),str(round(pi/2,rdigit)))}

device_numQubit = {'ibmqx2':5,
                    'ibmq_armonk':1,
                   'ibmq_athens':5,
                   'ibmq_vigo': 5,
                   'ibmq_ourense':5,
                   'ibmq_valencia':5,
                   'ibmq_santiago':5,
                   'ibmq_lima':5,
                   'ibmq_quito':5,
                   'ibmq_bogota': 5}
token0 = '316cca1ff42067b9901b6a9e90ca76496e76b4a73139954176c6fe31a45572fcd37089924cbe5c41edb96a9b4c0e92e44451ed350e2a7e45dc809d26cd752eed'
provider = IBMQ.enable_account(token0)
ddGateStr = {'X': 'u3(3.141592653589793,0,3.141592653589793)',
           'Y': 'u3(3.141592653589793,1.5707963267948966,1.5707963267948966)'}


def random_u3_gate():
    mat = operators.random.random_unitary(2)
    randu = unitary.UnitaryGate(mat)
    qstr = randu.qasm()
    theta, phi, lamb = [round(float(n),rdigit) for n in qstr.split(' ')[-2].split('(')[1].split(')')[0].split(',')]
    #TODO: a test function to make sure u_in and u_out are the same, can we use diamond norm?
    mat_out = library.U3Gate(theta,phi,lamb)
    phases = [mat.data.flatten()[i]/mat_out.to_matrix().flatten()[i] for i in range(4)]
    equal = all(abs(phases[0] - ele) < 0.1**rdigit for ele in phases)
    if equal:
        return ('u3',str(theta),str(phi),str(lamb))
    else:
        raise ValueError('U3 gate parameters does not match the input random unitary!')

prep_params = {"XplusState": ('u2',str(0),str(round(pi,rdigit))),
               'XminusState':('u2',str(round(pi,rdigit)),str(0)),
               'YplusState':('u2',str(round(pi/2,rdigit)),str(round(pi/2,rdigit))),
               'YminusState':('u2',str(round(-pi/2,rdigit)),str(round(-pi/2,rdigit))),
               'ZplusState':('id',''),
               'ZminusState':('u3',str(round(pi,rdigit)),str(-round(pi/2,rdigit)),str(round(pi/2,rdigit)))}
# construct 4 tetrahegron states
tetrahegrons_in_cartesian = {'tetra0': (np.sqrt(8/9),0,-1/3),
                             'tetra1': (-np.sqrt(2/9), np.sqrt(2/3), -1/3),
                             'tetra2': (-np.sqrt(2/9),-np.sqrt(2/3),-1/3),
                             'tetra3': (0,0,-1)}
def gate_params_for_unitary(mat):
    mat = Operator(np.array(mat))
    randu = unitary.UnitaryGate(mat)
    qstr = randu.qasm()
    angles_0 = [n for n in qstr.split(' ')[-2].split('(')[1].split(')')[0].split(',')]
    angles = [n.replace('pi', str(np.pi)) if 'pi' in n else n for n in angles_0]
    for i, n in enumerate(angles):
        if '*' in n and '/' not in n:
            numbers = angles[i].split('*')
            numbers = [float(f) for f in numbers]
            angles[i] = str(np.prod(numbers))
        elif '/' in n and '*' not in n:
            numbers = angles[i].split('/')
            numbers = [float(f) for f in numbers]
            angles[i] = str(numbers[0]/numbers[1])
        elif '/' in n and '*' in n:
            numbers = re.split('[*/]',n)
            numbers = [float(f) for f in numbers]
            m = n.index('*')
            d = n.index('/')
            if m<d:
                angles[i] = str(numbers[0] * numbers[1] / numbers[2])
            else:
                angles[i] = str(numbers[0] / numbers[1] * numbers[2])
    theta, phi, lamb = [round(myStr2float(n), rdigit) for n in angles]
    # TODO: a test function to make sure u_in and u_out are the same, can we use diamond norm?
    mat_out = library.U3Gate(theta, phi, lamb)
    phases = []
    for m,n in zip(mat.data.flatten(),mat_out.to_matrix().flatten()):
        if m <= 0.1**rdigit and n<=0.1**rdigit:
            phases.append(1)
        else:
            phases.append(m/n)
    equal = all(abs(phases[0] - ele) < 0.1 ** rdigit for ele in phases) # phases are the same for all elements
    if equal:
        return ('u3', str(theta), str(phi), str(lamb))
    else:
        raise ValueError('U3 gate parameters does not match the input random unitary!')

def gate_params_from_bloch_vector(bvector):
    """
    :param bvector: tuples (x,y,z) bloch vector in cartesian coordinates
    :return: gate name and gate params for state initialization
    """
    print(bvector)
    r, theta, phi = cartesian_to_polar(bvector)
    u00, u10 = np.cos(theta/2), np.sin(theta/2)* (np.cos(phi) + 1j * np.sin(phi))
    u01 = np.sqrt(u10 * np.conj(u10)/(u10 * np.conj(u10) + u00 * np.conj(u00)))
    u11 = - u01 * np.conj(u00) / np.conj(u10)
    ugate = np.array([[u00,u01],[u10,u11]])
    gstate = np.array([[1],[0]])
    psi = np.matmul(ugate, gstate)
    rho = np.matmul(psi,np.matrix.getH(psi))
    print('unitary')
    print([[u00,u01],[u10,u11]])
    print('bloch')
    print([np.trace(np.matmul(p,rho)) for p in [paulix,pauliy,pauliz]])
    return gate_params_for_unitary([[u00,u01],[u10,u11]])

tetrahegron_params = {key : gate_params_from_bloch_vector(item) for key,item in tetrahegrons_in_cartesian.items()}

prep_params.update(tetrahegron_params)

def construct_exp_dict(device,pstate,mbasis,circuitPath,**kwargs):
    """
    make the experiment configuration dict for free evolution
    """
    today = date.today()
    datestr = today.strftime("%Y%m%d")
    if 'spec_num' in kwargs:
        spec_num = kwargs.get('spec_num')
    else:
        spec_num = device_numQubit[device] - 1

    if 'runname' in kwargs:
        runname = kwargs.get('runname')
    elif 'measurement_config' in kwargs:
        runname = 'Meas%sFree'%(kwargs.get('measurement_config').capitalize())
    if 'mainq' in kwargs and isinstance(pstate,dict):
        runname = '_'.join([runname,'Q%d'%(int(kwargs.get('mainq'))),pstate['main'],'QS%d'%(spec_num),pstate['spec']])
    elif isinstance(pstate,str):
        runname = kwargs.get('runname')

    exp_dict = {"runname": runname, "date": datestr, "device": device, "spec_num":spec_num,
                 "circuitPath": circuitPath,'pstate': pstate, 'mbasis': mbasis,
                "filepath": circuitPath + device + '/' + datestr + '/' + runname + '/'}
    if 'mainq' in kwargs:
        exp_dict['mainq'] = int(kwargs.get('mainq'))
    if 'numIdGates' in kwargs:
        exp_dict['numIdGates'] = kwargs.get('numIdGates')
        exp_dict["filename"]= runname + '_' + datestr + '_' + device + '_' + 'numIdGates=' + str(
            exp_dict['numIdGates']) + '_' + exp_dict['mbasis'] + '.qasm'
    elif 'numDDrep' in kwargs:
        exp_dict['numDDrep'] = kwargs.get('numDDrep')
        exp_dict["filename"] = runname + '_' + datestr + '_' + device + '_' + 'numXYXYseq=' + str(
            exp_dict['numDDrep']) + '_' + mbasis + '.qasm'
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

def make_specq_list(n,reverse=False):
    if reverse:
        return list(range(4-n+1,5))
    else:
        return list(range(1, n+1))

def write_state_prep(f,exp_dict,**kwargs):
    """
    pstate is a string or a dictionary of strings
    """
    numQubits = device_numQubit[exp_dict['device']]
    if 'reverse' in kwargs:
        reverse = kwargs.get('reverse')
    else:
        reverse = False
    if 'mainq' in exp_dict:
        mainq = exp_dict['mainq']
    for i in range(numQubits):
        qubitStr = 'q[%d]' % (i)
        # get the gate type
        if i == mainq:
            gatename = exp_dict["pstate"]['main']
        elif i in make_specq_list(exp_dict['spec_num'],reverse=reverse):
            gatename = exp_dict["pstate"]['spec']
        else:
            gatename = 'ZplusState'
        if gatename[0] == 'randomState'[0]: # random unitary string
            rprep_params = kwargs.get('rprep_params')
            prepGateStr = 'u3(' + rprep_params[1] + ',' + rprep_params[2] + ',' + rprep_params[3] + ')'
        else:
            gate = prep_params[gatename][0] if gatename is not 'None' else None
            if gate is not None:
                if gate == 'u2':
                    prepGateStr = 'u2(' + prep_params[gatename][1] + ',' + prep_params[gatename][
                        2] + ')'
                elif gate == 'u3':
                    prepGateStr = 'u3(' + prep_params[gatename][1] + ',' + prep_params[gatename][
                        2] + ',' + prep_params[gatename][3] + ')'
                elif gate == 'id':
                    prepGateStr = 'id'
        f.write(prepGateStr + ' ' + qubitStr + ';\n')

def write_state_tomo_free_qasms(measure_config='all',**kwargs):
    """
    advance version of free evolution with option to prepare spectator qubits with user-specified initial states
    exp_dict[pstate] is a dictionary contains initial states for main qubit and spectator qubits
    """
    if "exp_dict" in kwargs:
        exp_dict = kwargs.get("exp_dict")
    else:
        print("input dictionary name error")
    Exp = TomoExperiment(exp_dict)
    if 'reverse' in kwargs:
        reverse = kwargs.get('reverse')
    else:
        reverse = False
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
    if exp_dict['pstate']['main'][0] == 'randomState'[0]:
        rprep_params = kwargs.get('rprep_params')
        write_state_prep(f, exp_dict,rprep_params=rprep_params,reverse=reverse)
    else:
        write_state_prep(f,exp_dict,reverse=reverse)
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
            if measure_config == 'mainq':
                if i == mainq:
                    qubitStr = 'q[%d]' %(i)
                    f.write(measGateStr + ' ' + qubitStr + ';\n')
            else:
                qubitStr = 'q[%d]' % (i)
                f.write(measGateStr + ' ' + qubitStr + ';\n')
    # measurement
    f.write(barrierStr)
    if measure_config == 'mainq':
        i = exp_dict['mainq']
        f.write('measure q[%d] -> c[%d];\n'%(i,i))
    else:
        for i in range(numQubits):
            f.write('measure q[%d] -> c[%d];\n' % (i, i))
    f.close()
    # transpilse and rewrite
    Exp.transpile_and_rewrite(Exp.filepath+Exp.filename)

class TomoExperiment:
    def __init__(self, exp_config):
        for key, value in exp_config.items():
            setattr(self, key, value)
    def get_backend(self):
        return getattr(provider.backends,self.device)
    def transpile_and_rewrite(self,qasm):
        """
        :param qasm: qasm filename
        :return:
        """
        backend = self.get_backend()
        f = QuantumCircuit().from_qasm_file(qasm)
        f = transpile(f, backend=backend, optimization_level=2)
        f2 = open(self.filepath + self.filename, "r+")
        f2.seek(0)
        f2.truncate()
        f2.write(f.qasm())
        f2.close()

def write_state_tomo_dd_on_spectator_qasms(measure_config='all',**kwargs):
    """
    DD on spec qubits on near qubits
    """
    if "exp_dict" in kwargs:
        exp_dict = kwargs.get("exp_dict")
    else:
        print("input dictionary name error")
    Exp = TomoExperiment(exp_dict)
    if 'reverse' in kwargs:
        reverse = kwargs.get('reverse')
    else:
        reverse = False
    if not os.path.exists(exp_dict["filepath"]):
        os.makedirs(exp_dict["filepath"])
    numQubits = device_numQubit[exp_dict['device']]
    mainq = int(exp_dict['mainq'])  # get the mainq index
    DDonSpecList = [mainq+1+2*i for i in range(2)]
    barrierStr = 'barrier %s;\n'%(write_qubit_index(numQubits))
    idStr = lambda q: 'id q[%s];\n' % (q)
    DDseq = ['X', 'Y', 'X', 'Y']
    f = open(exp_dict["filepath"] + exp_dict["filename"], "w")
    # header for 5 qubit
    f.write("OPENQASM 2.0;\ninclude\"qelib1.inc\";\nqreg q[%s];\ncreg c[%s];\n"%(numQubits,numQubits))
    # state preparation
    if exp_dict['pstate']['main'][0] == 'randomState'[0]:
        rprep_params = kwargs.get('rprep_params')
        write_state_prep(f, exp_dict, rprep_params=rprep_params, reverse=reverse)
    else:
        write_state_prep(f, exp_dict, reverse=reverse)
    # free evolution w. DD
    for d in range(exp_dict["numDDrep"]):
        for i in range(len(DDseq)): #XYXY sequence, or XY$
            f.write(barrierStr)
            for j in range(numQubits):
                if j in DDonSpecList:
                    f.write(ddGateStr[DDseq[i]] + ' q[%s];\n' % (j))
                else:
                    f.write(idStr(j))
    # pre-measurement rotation
    if meas_params[exp_dict["mbasis"]] is not None:
        f.write(barrierStr)
        measGateStr = 'u2('+meas_params[exp_dict["mbasis"]][0]+','+meas_params[exp_dict["mbasis"]][1]+')'
        for i in range(numQubits):
            if measure_config == 'mainq':
                if i == mainq:
                    qubitStr = 'q[%d]' % (i)
                    f.write(measGateStr + ' ' + qubitStr + ';\n')
            else:
                qubitStr = 'q[%d]' % (i)
                f.write(measGateStr + ' ' + qubitStr + ';\n')
    # measurement
    f.write(barrierStr)
    if measure_config == 'mainq':
        i = exp_dict['mainq']
        f.write('measure q[%d] -> c[%d];\n'%(i,i))
    else:
        for i in range(numQubits):
            f.write('measure q[%d] -> c[%d];\n'%(i,i))
    f.close()
    # transpilse and rewrite
    Exp.transpile_and_rewrite(Exp.filepath+Exp.filename)


circuitPath = r"../Circuits/"
device = "ibmq_bogota"
measurement_config = 'mainq'
# runname = 'MeasMainStateTomo_freeEvoLong'
measurement_basis = list(meas_params.keys())
# folder /device/date/state_tomography_freeEvo/
# pstates = ['XminusState', 'YminusState','YplusState']
#pstates = ['XplusState']
def write_free_and_dd():
    for pstate in pstates:
        runname = 'MeasAllstateTomo_freeEvo'
        runname = runname + '_' + pstate
        # measurement mitigation circuits
        """
        for p in ['ZplusState','ZminusState']:
            dict0 = construct_exp_dict(runname=runname,device=device,pstate=p,mbasis='obsZ',circuitPath=circuitPath)
            write_measurement_error_mitigation_qasm(p,'obsZ',exp_dict=dict0)
        """
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
        runname = 'MeasAllstateTomo_DD' + '_' + pstate
        # for p in ['ZplusState','ZminusState']:
        #     dict0 = construct_exp_dict(runname=runname,device=device,pstate=p,mbasis='obsZ',circuitPath=circuitPath)
        #     write_measurement_error_mitigation_qasm(p,'obsZ',exp_dict=dict0)
        for r in range(num_complete,num_repetition):
            numDDgates = sampling_rate * r
            for mbasis in measurement_basis:
                dict0 = construct_dd_exp_dict(runname=runname,device=device,pstate=pstate,mbasis=mbasis,circuitPath=circuitPath,numDDrep=numDDgates)
                write_state_tomo_dd_qasms(exp_dict=dict0)

def write_free_and_dd_w_spectators(qindex,pstate,spec_num=4,reverse=False, include_dd=False):
    """
    :param qindex:
    :param pstate:
    :param spec_num:
    :param reverse:
    :param include_dd:
    :return:
    """
    # # dense sampling rate
    num_repetition = 24
    num_complete = 0  # completed circuits
    sampling_rate = 144 # 48-72 for athens, 24 for ibmqx2, armonk,ibmqx2 vigo, 36 for valencia,72 for santiago
    rprep_params = random_u3_gate() if pstate['main'][0] == 'randomState'[0] else None
    """
    for r in range(num_complete, num_repetition):
        numIdGates = sampling_rate * r
        for mbasis in measurement_basis:
            dict0 = construct_exp_dict(runname='MeasMainqFreeLong',
                mainq=qindex, device=device, pstate=pstate, mbasis=mbasis,
                                       circuitPath=circuitPath, numIdGates=numIdGates,measurement_config=measurement_config,spec_num=spec_num)
            write_state_tomo_free_qasms(measure_config=measurement_config,exp_dict=dict0,rprep_params=rprep_params,reverse=reverse)
    # measurement mitigation circuits
    for p in ['ZplusState', 'ZminusState']:
        dict0 = construct_exp_dict(mainq=qindex,runname=dict0['runname'], device=device, pstate=p, mbasis='obsZ',
                                   circuitPath=circuitPath,spec_num=spec_num)
        write_measurement_error_mitigation_qasm(p, 'obsZ', exp_dict=dict0)
    """
    # DD circtuis
    if include_dd:
        for r in range(num_complete, num_repetition):
            numDDgates = sampling_rate * r
            for mbasis in measurement_basis:
                dict0 = construct_exp_dict(runname='MeasMainqFreeLongXY4',
                        mainq=qindex,device=device, pstate=pstate, mbasis=mbasis,
                                              circuitPath=circuitPath, numDDrep=int(numDDgates/4),measurement_config=measurement_config,spec_num=spec_num)
                write_state_tomo_dd_on_spectator_qasms(measure_config=measurement_config,exp_dict=dict0,rprep_params=rprep_params,reverse=reverse)
        for p in ['ZplusState', 'ZminusState']:
            dict0 = construct_exp_dict(mainq=qindex,runname=dict0['runname'], device=device, pstate=p, mbasis='obsZ',
                                       circuitPath=circuitPath,spec_num=spec_num)
            write_measurement_error_mitigation_qasm(p, 'obsZ', exp_dict=dict0)
pkeys5 = ['XplusState', 'XminusState', 'YplusState', 'YminusState',  'ZminusState']
pkeyr = 'randomState'
# a new config: 4 tetrahegron states and 1 random state
trkeys = list(tetrahegrons_in_cartesian.keys()) + [pkeyr]
nonMark_keys = ['XplusState','XminusState','tetra0','YplusState','YminusState']
spec_states = ['ZplusState','ZminusState','XplusState']
def main():
    for s in spec_states:
        for p in trkeys:
            pstates = {'main': p,
                       'spec': s}
            write_free_and_dd_w_spectators(qindex=0,pstate=pstates,spec_num=4, reverse=False, include_dd=True)


if __name__ == "__main__":
    main()
# for i in range(1):
#     pstates = {'main': pkeyr + str(i),
#                'spec': 'ZminusState'}
#     write_free_and_dd_w_spectators(qindex=3,pstate=pstates)

# check if a random qasm is valid
# qasm_name = dict0["filepath"] + dict0["filename"]
# qasm_obj = qasm.Qasm(qasm_name)
# qasm_obj.parse()
# #
# qasm_name2 = '../Circuits/ibmqx2/09112020/stateTomo_freeEvo_XplusState/' + 'stateTomo_freeEvo_XplusState_09112020_ibmqx2_numIdGates=48_XplusState_obsX.qasm'
# qasm_obj2 = qasm.Qasm(qasm_name2)
# qasm_obj2.parse()