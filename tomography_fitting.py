method = 'ibm'
import os
import pickle
import csv
import numpy as np
import matplotlib.pyplot as plt
from math import acos, atan
if method == 'rigetti':
    import qutip as qt
    from collections import namedtuple
    import glob
    from tomography.state_tomography import *
    from tomography.tomography import TomographySettings, DEFAULT_SOLVER_KWARGS
    from tomography.operator_utils import OperatorBasis


elif method == 'ibm':
    from IBM_QST_function import *
    from utility import hist2dict
    SX = np.array([[0, 1], [1, 0]])
    SY = np.array([[0, -1j], [1j, 0]])
    SZ = np.array([[1, 0], [0, -1]])

qubits = [0]  # the list of qubits to perform tomography on
nsamples = 8192
# an example datafile
datafile = r"stateTomo_freeEvo_11062020_ibmq_armonk_numIdGates=0_XplusState_obsX_5ee267acc926d60014453300.txt"
idGateTime = {'ibmqx2': 3.555555555555556e-02,
              'ibmq_armonk': 1.4222222222222224e-01}

def make_one_hist(datapath,qubiti,shots):
    """
    make histogram for the tomography experiment of one measurement observable
    :param datapath:
    :return:
    """
    file = open(datapath, "r")
    datalines = file.readlines()
    file.close()
    raw_data = []
    for line in datalines[2:]:
        bitstring = line.split(',')[0].split("\"")[1]
        count = int(line.split(',')[1].split("\n")[0])
        raw_data.append([bitstring, count])
    return count_outcome(raw_data,qubiti,shots)

def count_outcome(raw_data, qubiti,shots):
    """
    count measurement outcomes for q[i]
    :param qubit:
    :return: data[bitlist_to_int] = counts
    """
    data = [0, 0]
    for line in raw_data:
        if line[0][qubiti[0]] == '0':
            data[0] += line[1]
        elif line[0][qubiti[0]] == '1':
            data[1] += line[1]
        else:
            raise ValueError("bit string value not identified")
    if sum(data) != shots:
        raise ValueError('histgram counting error')
    return data



def locate_datafile(expobj,obs,numIdGates):
    dataname = expobj.runname + '_' + expobj.datestr + '_' + expobj.device + '_' + 'numIdGates=' + str(
        numIdGates) + '_' + expobj.pstate + '_' + obs
    data_fullpath = '/home/haimeng/LocalProjects/IBM-PMME/Data/raw/' + expobj.device + '/' + expobj.datestr + '/' + expobj.runname + '/' + expobj.runNo + '/' + dataname + '*.txt'
    datafile = glob.glob(data_fullpath)
    if len(datafile)>=1:
        return datafile[0]
    if len(datafile)<1:
        print(data_fullpath+' not found!')


def make_tomography_hist(expobj,numIdGates,qubiti,shots):
    """
    :param expobj:
    :param numIdGates:
    :param qubiti:
    :return: 3*2
    """
    mbasis = ['obsZ', 'obsX', 'obsY']
    hist = []
    for m in mbasis:
        datafile = locate_datafile(expobj,m,numIdGates)
        if datafile is not None:
            hist.append(make_one_hist(datafile,qubiti=qubiti,shots=shots))
    return np.array(hist)

pi = np.pi

rdigit = 5
# post-rotation gate
meas_params = {"obsZ": None,
               "obsX": ((round(pi, rdigit)), 0),
               "obsY": (0, str(round(3 * pi / 2, rdigit)))}
prep_params = {"XplusState": (str(round(pi,rdigit)),str(0))}

def U2Gate(angles):
    isqrt2 = 1 / np.sqrt(2)
    phi, lam = float(angles[0]), float(angles[1])
    u2matrix = np.array([
        [
            isqrt2,
            -np.exp(1j * lam) * isqrt2
        ],
        [
            np.exp(1j * phi) * isqrt2,
            np.exp(1j * (phi + lam)) * isqrt2
        ]
    ], dtype=complex)
    return qt.Qobj(u2matrix)



def to_density_matrix(state):
    return state * state.dag()

if method == 'rigetti':
    UNIT_TRACE = 'unit_trace'
    DEFAULT_STATE_TOMO_SETTINGS = TomographySettings(
        constraints={UNIT_TRACE},
        solver_kwargs=DEFAULT_SOLVER_KWARGS
    )

    # pauli
    QX = qt.sigmax()
    QY = qt.sigmay()
    QZ = qt.sigmaz()
    QI = qt.qeye(2)

    TOMOGRAPHY_GATES = [QI, U2Gate(meas_params["obsX"]), U2Gate(meas_params["obsY"])]

    DiagonalPOVM = namedtuple("DiagonalPOVM", ["pi_basis", "confusion_rate_matrix", "ops"])

    pauli_basis = OperatorBasis(
        [("I", QI / np.sqrt(2)), ("X", QX / np.sqrt(2)),
         ("Y", QY / np.sqrt(2)), ("Z", QZ / np.sqrt(2))])

    channel_ops = TOMOGRAPHY_GATES
    assignment_probs = np.array([[1, 0], [0, 1]])
    num_qubits = 1
    readout_povm = o_ut.make_diagonal_povm(o_ut.POVM_PI_BASIS ** num_qubits, assignment_probs)

def make_povm(confusion_rate_matrix=[[1, 0], [0, 1]]):
    confusion_rate_matrix = np.asarray(confusion_rate_matrix)
    gs = to_density_matrix(qt.basis(2, 0))
    es = to_density_matrix(qt.basis(2, 1))
    pi_basis = [gs, es]
    ops = [sum((pi_j * pjk for (pi_j, pjk) in zip(pi_basis, pjs)), 0)
           for pjs in confusion_rate_matrix]
    return DiagonalPOVM(pi_basis=pi_basis, confusion_rate_matrix=confusion_rate_matrix, ops=ops)

class StateTomographyFit:
    def __init__(self):
        self.idlist = None
        self.rhomat = None
        self.bloch = None
        self.pmat = None
        self.tlist = None
        self.runname = None
        self.config = None
        self.qindex = None

    def fit_state_from_run(self,datapath,exp0,qubiti):
        """
        fit qubit time evlotion tomography from raw data
        :param datapath: raw data path, str
        :param exp0: run info, ExpObject namedtuple
        :param qubiti: list of qubit index
        :return: None
        """
        self.qindex = qubiti
        self.config = exp0
        self.runname = '_'.join([exp0.runname, exp0.datestr , exp0.device, exp0.pstate, exp0.runNo])
        self.runname += '_Q%s'%(str(qubiti[0])) # single qubit tomography
        # change the information below
        datafiles_run = glob.glob(datapath+'*.txt')
        if '_' in exp0.device:
            strid = 5
        else:
            strid = 4
        idlist = sorted(list(set([int(f.split('/')[-1].split('_')[strid].split('=')[1]) for f in datafiles_run])))
        self.idlist = []
        self.rhomat = []
        self.pmat = []
        self.bloch = []
        def Qobj2BlochVec(qobj):
            vlist = [qobj * QX,qobj * QY,qobj * QZ]
            vector = np.array([a.tr() for a in vlist])
            return vector
        def rho2BlochVec(rho,polar=False):
            vlist = [np.matmul(rho,s) for s in [SX,SY,SZ]]
            vector = np.array([np.real(np.trace(a)) for a in vlist]) # Bloch vector has real values
            if polar:
                v = np.sqrt(np.sum(vector ** 2))
                theta = acos(vector[2]/v)
                phi = atan(vector[1]/vector[0]) if vector[0] != 0 else atan(np.inf)
                return np.array([v,theta,phi])
            else:
                return vector

        def rigetti_tomography(histograms):
            tomo_result = StateTomography.estimate_from_ssr(histograms, readout_povm, channel_ops,
                                                            settings=DEFAULT_STATE_TOMO_SETTINGS)
            return tomo_result
        def ibm_tomography(histograms):
            basisStrListIN = ['X', 'Y', 'Z']
            countsListIN = hist2dict(histograms,basisStrListIN)
            rho = IBM_QST_from_counts(countsListIN, basisStrListIN)
            return rho

        for i in range(len(idlist)):
            numIdGates = int(idlist[i])
            histograms = make_tomography_hist(exp0, numIdGates=numIdGates,qubiti=qubiti,shots=nsamples)
            if histograms.shape[0] == 3:
                self.idlist.append(numIdGates)
                self.pmat.append(histograms.flatten()/nsamples)
                if method == 'rigetti':
                    tomo_result = rigetti_tomography(histograms)
                    self.rhomat.append(qt.operator_to_vector(tomo_result.rho_est).data)
                    self.bloch.append(Qobj2BlochVec(tomo_result.rho_est))
                elif method == 'ibm':
                    tomo_result = ibm_tomography(histograms)
                    self.rhomat.append(tomo_result.transpose().flatten())
                    self.bloch.append(rho2BlochVec(tomo_result,polar=False))

            else:
                print("did not find all three config for tomography at id=%s"%(numIdGates))
        self.idlist = np.array(self.idlist)
        self.tlist = self.Gate2Time(gatetime=idGateTime[exp0.device])


        def Qobjlist2Nparray(Qlist):
            qmat = np.empty((4,len(Qlist)),dtype=complex)
            for i in range(len(Qlist)):
                qmat[:,i] = Qlist[i].todense().reshape((4,))
            return qmat

        if method == 'rigetti':
            self.rhomat = Qobjlist2Nparray(self.rhomat)
        elif method == 'ibm':
            self.rhomat = np.array(self.rhomat).transpose()
        self.bloch = np.array(self.bloch).transpose()
        self.pmat = np.array(self.pmat).transpose()

    def Gate2Time(self,gatetime):
        return self.idlist * gatetime

    def plot_bloch_vector(self,save=False,spacing=1,**kwargs):
        self.tlist = self.Gate2Time(gatetime=35)
        fig, ax = plt.subplots(ncols=1, nrows=1)
        fig.set_size_inches(7, 5)
        if 'polar' in kwargs:
            polar = kwargs.get('polar')
        if polar:
            labels = [r'$|v|$', r'$\theta$', r'$\phi$']
        else:
            labels = [r'$v_{x}$',r'$v_{y}$',r'$v_{z}$']
        for i in range(self.bloch.shape[0]):
            ax.plot(self.idlist[0::spacing], self.bloch[i,0::spacing],'-s',label=labels[i],color=mycolor[ckeys[i]],markersize=3,markeredgewidth=0.3, markeredgecolor=mycolor['black'])
        if 'idMax' in kwargs:
            idMax = kwargs.get('idMax')
        else:
            idMax = max(self.idlist)
        ax.set_xlim(-10,idMax)
        ax2 = ax.twiny()
        new_tick_locations = np.array([500*n for n in range(int(max(self.idlist)/500)+1)])
        def tick_function(X):
            V = X*idGateTime[self.config.device]
            return ["%.1f" % z for z in V]

        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(new_tick_locations)
        ax2.set_xticklabels(tick_function(new_tick_locations))
        ax2.set_xlabel(r"Evolution time ($\mu s$)")
        ax.set_xlabel(r"Number of identity gates")
        ax.legend(ncol=1, loc=4, prop={'size': 13}, frameon=True)
        plt.suptitle(r'Free evolution state tomography on %s Q%s $|-\rangle$ state'%(self.config.device,self.qindex[0]))
        plt.tight_layout(rect=[0,0,1,0.95])
        if save:
            if 'filepath' in kwargs:
                filepath = kwargs.get('filepath')
            if not os.path.exists(filepath):
                os.mkdir(filepath)
            filename = self.runname + '_%s'%(method) + '_blochVector'
            if polar:
                filename += '_polar'
            plt.savefig(filepath + filename + '.png')
            #pickle.dump(filepath, open(filepath + filename + '.pickle', "wb"))
            plt.show()
        else:
            plt.show()
        plt.close()


    def plot_denisty_matrix(self,save=False,**kwargs):
        self.tlist = self.Gate2Time(gatetime=35)
        fig, axs = plt.subplots(ncols=2, nrows=1)
        fig.set_size_inches(5 * 2, 5)
        labels = [r'$\rho_{00}$',r'$\rho_{10}$',r'$\rho_{01}$',r'$\rho_{11}$']
        titles = ['Real','Imaginary']
        for i in range(self.rhomat.shape[0]):
            axs[0].plot(self.idlist, np.real(self.rhomat)[i,:],'-o',label=labels[i],color=mycolor[ckeys[i]],markersize=5)
            axs[1].plot(self.idlist, np.imag(self.rhomat)[i,:],'-o',label=labels[i],color=mycolor[ckeys[i]],markersize=5)
        axs2 = [ax.twiny() for ax in axs]
        new_tick_locations = np.array([500*n for n in range(int(max(self.idlist)/500)+1)])
        def tick_function(X):
            V = X*idGateTime[self.config.device]
            return ["%.1f" % z for z in V]
        for i in range(len(axs2)):
            axs2[i].set_xlim(axs[i].get_xlim())
            axs2[i].set_xticks(new_tick_locations)
            axs2[i].set_xticklabels(tick_function(new_tick_locations))
            axs2[i].set_xlabel(r"Evolution time ($\mu s$)")
            axs[i].set_xlabel(r"Number of identity gates")
            axs[i].legend(ncol=2, loc=4, prop={'size': 13}, frameon=True)
            axs[i].set_title(titles[i])
        plt.suptitle(r'Free evolution state tomography on %s Q%s $|-\rangle$ state'%(self.config.device,self.qindex[0]),y=1.1)
        plt.tight_layout(rect=[0,0,1,0.95])
        if save:
            if 'filepath' in kwargs:
                filepath = kwargs.get('filepath')
            if not os.path.exists(filepath):
                os.mkdir(filepath)
            filename = self.runname
            plt.savefig(filepath + filename + '.png')
            #pickle.dump(filepath, open(filepath + filename + '.pickle', "wb"))
            plt.show()
        else:
            plt.show()
        plt.close()

    def saveBloch2csv(self,filepath,**kwargs):
        data = zip(self.idlist,self.tlist,np.real(self.bloch[0,:]),np.real(self.bloch[1,:]),np.real(self.bloch[2,:]))
        filename = self.runname +'_%s'%(method) +'_blochVector'
        if 'polar' in kwargs:
            polar = kwargs.get('polar')
        else:
            polar = False
        if polar:
            header = ['id', 'time', '|v|', 'theta', 'phi']
            filename += '_polar'
        else:
            header = ['id', 'time', 'vx', 'vy', 'vz']
        with open(filepath + filename + '.csv', 'w', newline='') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(header)
            writer.writerows(data)

    def saveRho2csv(self,filepath):
        data = zip(self.idlist,self.tlist,np.real(self.rhomat[0,:]),np.imag(self.rhomat[0,:]),np.real(self.rhomat[1,:]),np.imag(self.rhomat[1,:]),\
                   np.real(self.rhomat[2,:]),np.imag(self.rhomat[2,:]),np.real(self.rhomat[3,:]),np.imag(self.rhomat[3,:]))
        filename = self.runname
        header = ['id','time','rho00(real)','rho00(imag)','rho10(real)','rho10(imag)','rho01(real)','rho01(imag)','rho11(real)','rho11(imag)']
        with open(filepath + filename+'_%s'%(method)  +'_densityMatrix.csv', 'w', newline='') as file:
            writer = csv.writer(file, delimiter=',')
            writer.writerow(header)
            writer.writerows(data)


num_qubits = len(qubits)
dimension = 2 ** num_qubits

ExpObject = namedtuple('ExpObject', ['runname', 'datestr', 'device', 'pstate', 'runNo'])

mycolor={ 'reddish purple':'#cc79a7',\
         'blush green': '#009e73', \
       'sky blue':'#56b4e9',\
        'orange':'#e69f00', \
        'vermillion':'#d55e00',\
       'blue':'#0072b2',\
        'yellow':'#f0e442',\
        'black':'#000000'}
ckeys = list(mycolor.keys())


# edit from below
# read date from
datapath = r'/home/haimeng/LocalProjects/IBM-PMME/Data/raw/'
# raw data info
exp0 = ExpObject(runname='stateTomo_freeEvo', datestr='09052020', device='ibmqx2', pstate='XplusState',
                 runNo='run-sum')
datapath = r"/home/haimeng/LocalProjects/IBM-PMME/Data/raw/" + exp0.device + "/" + exp0.datestr + "/" + exp0.runname + "/" + exp0.runNo +"/"
# store date to
filepath = r"/home/haimeng/LocalProjects/IBM-PMME/Analysis/" + exp0.device + "/" + exp0.datestr + '/'
tomo_fit = StateTomographyFit()
for i in range(5):
    tomo_fit.fit_state_from_run(datapath,exp0,qubiti=[i])

    # change where x-axis stops using keyword argument 'idMax';
    # change x-axis resolution using argument 'spacing', sampling frequency equals 1 sample per spacing*4 Id gates
    tomo_fit.plot_bloch_vector(save=True,filepath=filepath,idMax=1200,spacing=1,polar=False)
    #tomo_fit.plot_denisty_matrix(save=False,filepath=filepath,spacing=2)
    tomo_fit.saveBloch2csv(filepath,polar=False)
    # tomo_fit.save2csv(filepath)
    # tomo_fit.saveRho2csv(filepath)
    # import pickle
    # filepath = r"/home/haimeng/LocalProjects/IBM-PMME/Analysis/"
    # filename = exp0.runname + '_' + exp0.datestr +'_'+exp0.device +'_'+exp0.pstate +'_run1.pickle'
    # pickle.dump(data, open( filepath+filename, "wb" ))