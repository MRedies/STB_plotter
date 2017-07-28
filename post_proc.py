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

def plot_square_lattice(root, linesty='', 
                        figsize=(8,6), axis=None, linelabel=None,
                        n_bands=8):
    if axis == None:
	fig, axis = plt.subplots(figsize=figsize)
    else: 
	fig = None
    #data = np.load(filename)
    E = np.load(root + 'band_E.npy')
    E = np.transpose(E)
    n_ktps = (E.shape[0]+2)/3
    ticks = calc_tick_pos(3, n_ktps)
    #print( ticks)
    label = ["$\Gamma$", "X", "M", "$\Gamma$"]

    try:
        E_fermi = np.load(root + 'E_fermi.npy')[0]
    except:
        E_fermi = None
    
    #print(E.shape)
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
    #axis.set_ylabel(get_unit(data, 'energy'), fontsize=15)
    axis.set_ylim([np.min(E), np.max(E)])
    axis.xaxis.grid(True)

    return fig, axis

def plot_hex_lattice(root, linesty='', ylim=None, 
                        figsize=(8,6), axis=None, linelabel=None):
    if axis == None:
	fig, axis = plt.subplots(figsize=figsize)
    else: 
	fig = None
    E = np.load(root + 'band_E.npy')
    E = np.transpose(E)
    
    if ylim == None:
        ylim = np.array([np.min(E), np.max(E)])
    sel = np.logical_and(np.any(E >= ylim[0], axis=0) 
                       , np.any(E <= ylim[1], axis=0))
    E = E[:,sel]
    
    n_ktps = (E.shape[0]+2)/3
    ticks = calc_tick_pos(3, n_ktps)
    #print( ticks)
    label = ["$\Gamma$", "M", "K", "$\Gamma$"]

    try:
        E_fermi = np.load(root + 'E_fermi.npy')[0]
    except:
        E_fermi = None
    
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
    #axis.set_ylabel(get_unit(data, 'energy'), fontsize=15)
    axis.set_ylim([np.min(E), np.max(E)])
    axis.xaxis.grid(True)

    return fig, axis

def inte_dos(dos, E):
    summe = 0.0
    dE = E[1] - E[0]
    for i in range(len(E)-1):
        summe += 0.5 * dE * (dos[i+1]  + dos[i])
    return summe

def check_dos(root):
    get = lambda name: np.load(root + name + '.npy')
    n_atm = 1.0 * get('DOS_partial').shape[0]
    
    n_dos = inte_dos(get('DOS'), get('DOS_E'))
    print("n_dos = {}".format(n_dos))
    if(np.abs((n_atm - n_dos)/(n_atm)) > 0.05):
        print "False gen_dos"
    else:
        pass
        print "Good"
    
    n_up = inte_dos(get('DOS_up'), get('DOS_E'))
    print("n_up = {}".format(n_up))
    if(np.abs((0.5*n_atm - n_up)/(0.5*n_atm)) > 0.05):
        print "False up_dos"
    else:
        pass
        print "Good"
    
    n_down = inte_dos(get('DOS_down'), get('DOS_E'))
    print("n_down = {}".format(n_down))
    if(np.abs((0.5*n_atm - n_down)/(0.5*n_atm)) > 0.05):
        print "False down_dos"
    else:
        pass
        print "Good"
        
    PDOS = get('DOS_partial')
    for i in range(get('DOS_partial').shape[0]):
        n = inte_dos(PDOS[i,:], get('DOS_E'))
        print("N = {}".format(n))
        if(np.abs(n - 1.0) > 0.05):
            print "DOS {} false".format(i)
        else:
            pass
            print "Good"


def get_mag_overall(folder, Ef):
    DOS_E    = np.load(folder + "DOS_E.npy")
    DOS_up   = np.load(folder + "DOS_up.npy")
    DOS_down = np.load(folder + "DOS_down.npy")

    sel  = DOS_E <= Ef
    up   = inte_dos(DOS_up[sel], DOS_E[sel])
    down = inte_dos(DOS_down[sel], DOS_E[sel])

    return up - down

def y_cut_idx(folder):
    x = np.load(folder + "pos_x.npy")
    y = np.load(folder + "pos_y.npy")
    return np.argwhere(np.logical_and((np.abs(x) < 1e-6), y < 0.0)) 


def plot_mag_cut(folder, figsize=(8,6), axis=None):
    cut = y_cut_idx(folder)
    x   = np.load(folder + "pos_x.npy")[cut]
    y   = np.load(folder + "pos_y.npy")[cut]
    theta = np.load(folder + "m_theta.npy")[cut]
    phi = np.load(folder + "m_phi.npy")[cut]

    r = 1.0
    m_y = np.sin(theta) * np.sin(phi)
    m_z = np.cos(theta)

    plt.quiver(-y, x, m_y, m_z, pivot="mid", scale=12)
    plt.plot(-y, x,"r.")
    plt.axes().set_aspect('equal', 'datalim')
    plt.show()


