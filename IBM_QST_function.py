
"""
Quantum State Tomography from IBM, using counts and measurement basis
Author: Matthew Kowalsky
Date: 6/30/20

Input: List of count dictionaries: [{'000': 435, '110': 217, ... },{'000': 321, ...},..]
       List of basis strings, each matched to a count dictionary ['XXX','XYX',...]
       Both lists should be length 3**n, and measurements in Pauli basis
       
Output: rho as a numpy.ndarray, where rho is an approximate density matrix of the system

"""
from qiskit import QuantumCircuit
import qiskit.result.result
import qiskit.result.models
import qiskit.validation.base
from qiskit.ignis.verification.tomography import StateTomographyFitter
from utility import *

global AllStates
AllStates = []

def returnAllKLengthStates(setL, k): 
    n = len(setL) 
    returnAllKLengthStatesRec(setL, "", n, k) 

'''Recursive method to return all strings of length k from set'''
def returnAllKLengthStatesRec(setL, prefix, n, k): 
    if (k == 0) : 
        AllStates.append(prefix)
        return
    for i in range(n): 
        newPrefix = prefix + setL[i] 
        returnAllKLengthStatesRec(setL, newPrefix, n, k - 1)
        

def IBM_QST_from_counts(countsListIN,basisStrListIN):
    num_qubits=len(basisStrListIN[0])
    num_shots= sum([val for val in countsListIN[0].values()])
    
    #Convert data into the form IBM's StateTomographyFitter() takes.
    
    #convert the list of strings like 'XXX' to a list of tuples like (X,X,X)
    #tupleListIN=[(elem[0],elem[1],elem[2]) for elem in basisStrListIN]
    tupleListIN = [(elem[0],) for elem in basisStrListIN]
    
    #Create list of all possible readout states as binary strings ex. ['001','000',....]
    returnAllKLengthStates(['0','1'],num_qubits) #Makes AllStates list complete.
    
    #convert the list of count dictionaries from having binary keys to hexadecimal keys
    hexaDecimal_countsList=[]
    for i in range(len(countsListIN)):
        hex_countsDict={}
        for key in AllStates:
            hex_key=str(hex(int(key, 2))) #convert key from binary to hex
            #print(key," ",hex_key) #A check
            hex_countsDict[hex_key]=countsListIN[i][key] #fill up new copy of dict with hex_keys
            #hex_countsDict={'0x0': 456, ... '0x7': 987}
            
            Obj_hex_countsDict=qiskit.validation.base.Obj(**hex_countsDict)#New IBM object, still only has hex_dict info
            data=qiskit.result.models.ExperimentResultData(**{'counts':Obj_hex_countsDict})#Higher level IBM object, still only has hex_dict info
            
        headerDictObj=qiskit.validation.base.Obj(**{'name': str(tupleListIN[i]), 'meas_level': 2, 'memory':False, 'n_qubits' : num_qubits, 'memory_slots': num_qubits })
        hexaDecimal_countsList.append(qiskit.result.models.ExperimentResult(**{'shots':num_shots, 'success':True, 'data': data, 'header': headerDictObj}))    
    #print(hexaDecimal_countsList)
    
    #Create the fake Result object and the fake List of QuantumCircuit objects
    Result_init_dict={'backend_name': 'bogus', 'job_id': 'bogus', 'success': True, 'results': hexaDecimal_countsList , 'qobj_id': 'bogus', 'backend_version': '2.0.5'}
    FakeResults=qiskit.result.result.Result(**Result_init_dict)
    FakeCircuits=[QuantumCircuit(num_qubits,num_qubits,name=str(tupleListIN[i])) for i in range(len(tupleListIN))]
    
    #Get the density matrix from StateTomographyFitter()
    tomoFAKED=StateTomographyFitter(FakeResults,FakeCircuits)
    rho = tomoFAKED.fit('auto') #methods are 'cvx','lstsq','auto'
    return rho #rho is a numpy.ndarray

# basisStrListIN = ['X','Y','Z']
# exp0 = ExpObject(runname='stateTomo_freeEvo', datestr='11062020', device='ibmqx2', pstate='XplusState',
#                  runNo='run1')
#
# countsListIN = returnTomoCountDict(exp0,0)
# rho = IBM_QST_from_counts(countsListIN,basisStrListIN)
