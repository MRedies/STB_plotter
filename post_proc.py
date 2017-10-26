# import matplotlib as mpl
# mpl.use('Qt4Agg')
from mpl_toolkits.mplot3d import Axes3D

import configparser
#import StringIO
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.tri as tri

def calc_kpts(n_sec, pts_per_sec):
    return n_sec * (pts_per_sec - 1) + 1

def calc_tick_pos(n_sec, pts_per_sec):
    ticks = np.zeros(n_sec+1)
    for i in range(n_sec + 1):
        ticks[i] = (pts_per_sec - 1) * i
    return ticks


def get_unit(folder, name):
    try:
        with open(folder + "setup.cfg") as f:
            lines = f.readlines()
            for l in lines:
                if name in l:
                    l = l.replace('"', "'")
                    return l.split("'")[1]
    except Exception as e:
        print(e)
        print("cant find unit 1")
        return " "
    print("cant find unit 2")
    return " "


def plot_square_lattice(root, linesty='', 
        figsize=(8,6), axis=None, linelabel=None):
    if axis is None:
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
    k = np.arange(E.shape[0])
    if linelabel is None:
        axis.plot(k,E, linesty)
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi)
    else:
        axis.plot(k,E, linesty, label=linelabel)
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi)

    axis.set_xticks(ticks)
    axis.set_xticklabels(label)
    axis.set_ylabel(get_unit(root, 'energy'), fontsize=25)
    axis.set_ylim([np.min(E), np.max(E)])
    axis.xaxis.grid(True)

    return fig, axis

def plot_hex_lattice(root, linesty='', ylim=None, linewidth=1, fontsize=16,
        figsize=(8,6), axis=None, linelabel=None):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)
    else: 
        fig = None
    E = np.load(root + 'band_E.npy')
    E = np.transpose(E)

    if ylim is None:
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

    if linelabel is None:
        axis.plot(k,E, linesty, linewidth=linewidth)
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi)
    else:
        axis.plot(k,E, linesty, label=linelabel, linewidth=linewidth)
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi)
    axis.tick_params(labelsize=fontsize)
    axis.set_xticks(ticks)
    axis.set_xticklabels(label)
    axis.set_ylabel(get_unit(root, 'energy'), fontsize=25)
    axis.set_ylim(ylim)
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
    all_good = True
    n_dos = inte_dos(get('DOS'), get('DOS_E'))
    #print("n_dos = {}".format(n_dos))
    if(np.abs((n_atm - n_dos)/(n_atm)) > 0.05):
        print("False gen_dos")
        all_good = False
    else:
        pass
    #print "Good"

    n_up = inte_dos(get('DOS_up'), get('DOS_E'))
    #print("n_up = {}".format(n_up))
    if(np.abs((0.5*n_atm - n_up)/(0.5*n_atm)) > 0.05):
        print("False up_dos")
        all_good = False
    else:
        pass
    #print "Good"

    n_down = inte_dos(get('DOS_down'), get('DOS_E'))
    #print("n_down = {}".format(n_down))
    if(np.abs((0.5*n_atm - n_down)/(0.5*n_atm)) > 0.05):
        print("False down_dos")
        all_good = False
    else:
        pass
    #print "Good"

    PDOS = get('DOS_partial')
    for i in range(get('DOS_partial').shape[0]):
        n = inte_dos(PDOS[i,:], get('DOS_E'))
        #print("N = {}".format(n))
        if(np.abs(n - 1.0) > 0.05):
            print("DOS {} false".format(i))
            all_good = False
        else:
            pass
    print("All good = {}".format(all_good))

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
    return np.argwhere(np.abs(x) < 1e-6)
#return np.argwhere(np.logical_and((np.abs(x) < 1e-6), y < 0.0)) 


def plot_mag_cut(folder, figsize=(8,6), axis=None, fontsize=16):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)
    else: 
        fig = None
    cut = y_cut_idx(folder)
    x   = np.load(folder + "pos_x.npy")[cut]
    y   = np.load(folder + "pos_y.npy")[cut]
    theta = np.load(folder + "m_theta.npy")[cut]
    phi = np.load(folder + "m_phi.npy")[cut]


    r = 1.0
    m_y = np.sin(theta) * np.sin(phi)
    m_z = np.cos(theta)

    axis.quiver(-y, x, m_y, m_z, pivot="mid", scale=12)
    axis.plot(-y, x,"r.")
    axis.set_aspect('equal', 'datalim')
    axis.tick_params(labelsize=fontsize)


def calc_local_dos(folder, E_window):
    pdos = np.load(folder + "DOS_partial.npy")
    Edos = np.load(folder + "DOS_E.npy")
    cut = np.argwhere(np.logical_and(Edos >= E_window[0], Edos <= E_window[1]))
    try:
        l_cut = np.min(cut)
        u_cut = np.max(cut)+1
    except:
        print("No data in Energy window")
        return None

    Edos = Edos[l_cut:u_cut]


    n_band = pdos.shape[0]
    n_loc = n_band /2
    loc_dos = pdos[:n_loc,l_cut:u_cut] + pdos[n_loc:,l_cut:u_cut]

    integr = np.zeros(loc_dos.shape[0])
    for i in range(loc_dos.shape[0]):
        integr[i] = inte_dos(loc_dos[i,:], Edos)
    return integr

def plot_loc_dos(folder, E_window, figsize=(8,8), axis=None, fig=None, n_lvls=20,
        lvls=None, cmap=plt.cm.viridis):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)

    loc_dos = calc_local_dos(folder, E_window)
    x = np.load(folder + "pos_x.npy")
    y = np.load(folder + "pos_y.npy")
    if(lvls is None):
        lvls = np.linspace(np.min(loc_dos), np.max(loc_dos), n_lvls)

    triang = tri.Triangulation(x, y)
    cbar = axis.tricontourf(triang, loc_dos,levels=lvls, cmap=cmap)
    axis.set_aspect('equal')

    cut = y_cut_idx(folder)
    axis.plot(x[cut], y[cut], 'g', linewidth=1)
    axis.plot(x,y, 'r.', markersize=2)

    fig.colorbar(cbar, ax=axis)

def plot_loc_dos_surf(folder, E_window, figsize=(8,8), axis=None, 
        fig=None, n_lvls=20, lvls=None, cmap=plt.cm.coolwarm):
    if axis is None:
        fig = plt.figure()
        axis = fig.gca(projection='3d')

    loc_dos = calc_local_dos(folder, E_window)
    x = np.load(folder + "pos_x.npy")
    y = np.load(folder + "pos_y.npy")
    if(lvls is None):
        lvls = np.linspace(np.min(loc_dos), np.max(loc_dos), n_lvls)

    axis.plot_trisurf(x,y,loc_dos, cmap=cmap,
            linewidth=0, antialiased=False)

    def plot_hall(folder, figsize=(8,8), axis=None, xlim=None, ylim=None, fontsize=16, linesty='', linewidth=1, label="$\sigma_{xy}$", color=None):
        if axis is None:
            fig, axis = plt.subplots(figsize=figsize)

    E = np.load(folder + "hall_cond_E.npy")
    c = np.load(folder + "hall_cond.npy")

    axis.plot(E,c, linesty, linewidth=linewidth, label=label, color=color)
    axis.set_ylabel(r"$\sigma_{xy}$", fontsize=fontsize)
    axis.set_xlabel(r"eV", fontsize=fontsize)
    axis.set_title(r"Hall conductance", fontsize=fontsize)

    if(not(xlim is None)):
        axis.set_xlim(xlim)
    if(not(ylim is None)):
        axis.set_ylim(ylim)

    axis.tick_params(labelsize=fontsize)
    axis.grid(True)


def plot_orbmag(folder, figsize=(8,8), axis=None, xlim=None, ylim=None, labels=["L", "IC", "M"],
        fontsize=16, linesty='', linewidth=1, which="ALL", color=None):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)

    E = np.load(folder + "orbmag_E.npy")
    if(xlim is None):
        sel = np.arange(E.shape[0])
    else:
        sel = np.argwhere(np.logical_and(E>=xlim[0], E<=xlim[1]))

    if(which.upper() == "M" or which.upper()=="ALL"):
        om = np.load(folder + "orbmag.npy")
        axis.plot(E[sel],om[sel], linesty, linewidth=linewidth, label=labels[2], color=color)
    if(which.upper() == "L" or which.upper()=="ALL"):
        L = np.load(folder + "orbmag_L.npy")
        axis.plot(E[sel],L[sel], linesty, linewidth=linewidth, label=labels[0], color=color)
    if(which.upper() == "IC" or which.upper()=="ALL"):
        IC = np.load(folder + "orbmag_IC.npy")
        axis.plot(E[sel],IC[sel], linesty, linewidth=linewidth, label=labels[1], color=color)

    axis.set_ylabel(r"M/$\mu_b$", fontsize=fontsize)
    axis.set_xlabel(r"eV", fontsize=fontsize)
    axis.set_title(r"Orbital magnetism", fontsize=fontsize)
    axis.tick_params(labelsize=fontsize)

    if(not(xlim is None)):
        axis.set_xlim(xlim)
    if(not(ylim is None)):
        axis.set_ylim(ylim)

def N_sk(fol, dims=None):
    if(dims is None ):
        s = np.load(fol + "pos_x.npy").shape[0]
        sqr = int(np.sqrt(s))
        if(sqr*sqr == s):
            dims = (sqr, sqr)
        else:
            raise ValueError("need to specifiy non-square dims")
                
    r = np.zeros((dims[0], dims[1],3))
    r[:,:,0]     = np.load(fol + "pos_x.npy").reshape(dims, order="F")
    r[:,:,1]     = np.load(fol + "pos_y.npy").reshape(dims, order="F")
    r[:,:,2]     = np.load(fol + "pos_z.npy").reshape(dims, order="F")

    phi   = np.load(fol + "m_phi.npy").reshape(dims, order="F")
    theta = np.load(fol + "m_theta.npy").reshape(dims, order="F")

    m = np.zeros((dims[0], dims[1],3))
    m[:,:,0] = np.sin(theta) * np.cos(phi)
    m[:,:,1] = np.sin(theta) * np.sin(phi)
    m[:,:,2] = np.cos(theta)

    g_x = np.gradient(m, 1.0, axis=0)
    g_y = np.gradient(m, 1.0, axis=1)
    cro = np.cross(g_x, g_y)

    summe = 0.0
    for i in range(40):
        for j in range(40):
            summe += np.inner(m[i,j,:], cro[i,j,:])
    summe *= 1.0/(4.0*np.pi)
    return summe
