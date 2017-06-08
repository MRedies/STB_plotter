# import matplotlib as mpl
# mpl.use('Qt4Agg')
import ConfigParser
import StringIO
import numpy as np
import matplotlib.pyplot as plt
import re

def calc_kpts(n_sec, pts_per_sec):
    return n_sec * (pts_per_sec - 1) + 1

def calc_tick_pos(n_sec, pts_per_sec):
    ticks = np.zeros(n_sec+1)
    for i in range(n_sec + 1):
        ticks[i] = (pts_per_sec - 1) * i
    return ticks


def get_unit(data, name):
    p = re.compile("\s*?energy\s*?=\s*?'(.*?)'")
    cfg = data['input.cfg']
    for line in cfg.split("\n"):
        if p.match(line) != None:
            return p.search(line).group(0).split("'")[1]

def plot_square_lattice(data, linesty='', 
                        fig_sz=(8,6), axis=None, linelabel=None,
                        n_bands=8):
    if axis == None:
	fig, axis = plt.subplots(figsize=fig_sz)
    else: 
	fig = None
    #data = np.load(filename)
    n_ktps = data['band_num_kpts'][0]
    ticks = calc_tick_pos(3, n_ktps)
    label = ["$\Gamma$", "X", "M", "$\Gamma$"]
    E = data['band_E']
    E = np.transpose(E)

    if("E_fermi" in data.files):
        E_fermi = data['E_fermi'][0]
    else: 
        E_fermi = None
    
    print(E.shape)
    E = E[:,:n_bands]
    k = np.arange(E.shape[0])
    if linelabel == None:
	axis.plot(k,E, linesty)
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi)
    else:
	axis.plot(k,E, linesty, label=linelabel)
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi)

    axis.set_xticks(ticks)
    axis.set_xticklabels(label)
    axis.set_ylabel(get_unit(data, 'energy'), fontsize=15)
    axis.set_ylim([np.min(E), np.max(E)])
    axis.xaxis.grid(True)

    return fig, axis

def inte_dos(dos, E):
    summe = 0.0
    dE = E[1] - E[0]
    for i in range(len(E)-1):
        summe += 0.5 * dE * (dos[i+1]  + dos[i])
    return summe

def check_dos(data):
    n_atm = 1.0 * data['DOS_partial'].shape[0]
    
    n_dos = inte_dos(data['DOS'], data['DOS_E'])
    #print("n_dos = {}".format(n_dos))
    if(np.abs((n_atm - n_dos)/(n_atm)) > 0.05):
        print "False gen_dos"
    
    n_up = inte_dos(data['DOS_up'], data['DOS_E'])
    #print("n_up = {}".format(n_up))
    if(np.abs((0.5*n_atm - n_up)/(0.5*n_atm)) > 0.05):
        print "False up_dos"
    
    n_down = inte_dos(data['DOS_down'], data['DOS_E'])
    #print("n_down = {}".format(n_down))
    if(np.abs((0.5*n_atm - n_down)/(0.5*n_atm)) > 0.05):
        print "False down_dos"
        
    PDOS = data['DOS_partial']
    for i in range(data['DOS_partial'].shape[0]):
        n = inte_dos(PDOS[i,:], data['DOS_E'])
        #print("N = {}".format(n))
        if(np.abs(n - 1.0) > 0.05):
            print "DOS {} false".format(i)
