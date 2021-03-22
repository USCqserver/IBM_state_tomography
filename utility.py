"""
Haimeng Zhang
to work with Matt's IBM tomography fitter
"""
from collections import namedtuple
import glob
import numpy as np
from math import atan2, acos
ExpObject = namedtuple('ExpObject', ['runname', 'datestr', 'device', 'pstate', 'runNo'])

def locate_datafile(expobj,obs,numIdGates):
    dataname = expobj.runname + '_' + expobj.datestr + '_' + expobj.device + '_' + 'numIdGates=' + str(
        numIdGates) + '_' + expobj.pstate + '_' + obs
    data_fullpath = '/home/haimeng/LocalProjects/IBM-PMME/Data/raw/' + expobj.device + '/' + expobj.datestr + '/' + expobj.runname + '/' + expobj.runNo + '/' + dataname + '*.txt'
    datafile = glob.glob(data_fullpath)
    if len(datafile)>=1:
        return datafile[0]
    if len(datafile)<1:
        print(data_fullpath+' not found!')

def countWithBinaryKey(datapath):
    """
    make histogram for the tomography experiment of one measurement observable
    :param datapath:
    :return:
    """
    file = open(datapath, "r")
    datalines = file.readlines()
    file.close()
    count_dict = {}
    for line in datalines[2:]:
        bitstring = line.split(',')[0].split("\"")[1]
        count = int(line.split(',')[1].split("\n")[0])
        count_dict[bitstring] = count
    return count_dict

def returnTomoCountDict(expobj,numIdGates):
    mbasis = ['obsX', 'obsY', 'obsZ']
    hist = []
    for m in mbasis:
        datafile = locate_datafile(expobj, m, numIdGates)
        if datafile is not None:
            hist.append(countWithBinaryKey(datafile))
    return hist

def hist2dict(histograms,basisStrListIN):
    """
    make sure the histogram is X,Y,Z ordered
    """
    sort_hist = {basisStrListIN[0]: histograms[1] , basisStrListIN[1]: histograms[2], basisStrListIN[2]: histograms[0]}
    hist_dict = [{'1': int(sort_hist[key][1]),'0':int(sort_hist[key][0])} for key in sort_hist.keys()]
    return hist_dict

def blochInPolarCoordinate(vector):
    v = np.sqrt(np.sum(vector ** 2))
    theta = acos(vector[2] / v)/np.pi
    phi = atan2(vector[1] , vector[0])/np.pi
    return np.array([v, theta, phi])

def cartesian_to_polar(bvector):
    """
    :param bvector: tuples (x,y,z)
    :return: polar coordinates (r, theta, phi) in radians
    """
    x,y,z = bvector
    r = np.sqrt(x**2 + y**2 + z**2)
    if x!=0:
        phi = np.arctan(y/x)
    elif x==0 and y!=0:
        phi = np.arctan(np.inf)
    elif x==0 and y==0:
        phi=0
    return (r,np.arccos(z/r),phi)


def myStr2float(string):
    if '/' in string:
        a,b = string.split('/')
        return float(a)/float(b)
    else:
        return float(string)