# import matplotlib as mpl
# mpl.use('Qt4Agg')
from mpl_toolkits.mplot3d import Axes3D

import configparser
#import StringIO
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.tri as tri
import matplotlib.image as mpimg
from glob import glob
from itertools import combinations

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
        figsize=(8,6), axis=None, linelabel=None, color=None):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)
    else:
        fig = None
    #data = np.load(filename)
    E = np.load(root + 'band_E.npy')
    E = np.transpose(E)
    n_ktps = (E.shape[0]+2)/3
    ticks = calc_tick_pos(3, n_ktps)
    label = ["$\Gamma$", "X", "M", "$\Gamma$"]

    try:
        E_fermi = np.load(root + 'E_fermi.npy')[0]
    except:
        E_fermi = None

    k = np.arange(E.shape[0])
    if linelabel is None:
        axis.plot(k,E, linesty, color=color)
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi, color=color)
    else:
        axis.plot(k,E, linesty, label=linelabel, color=color)
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi, color=color)

    axis.set_xticks(ticks)
    axis.set_xticklabels(label)
    axis.set_ylabel(get_unit(root, 'energy'), fontsize=25)
    dE = np.max([np.max(E) - np.min(E), 0.01])
    axis.set_ylim([np.min(E)-0.03*dE, np.max(E)+0.03*dE])
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
        axis.plot(k,E[:,0], linesty, label=linelabel, linewidth=linewidth)
        axis.plot(k,E[:,0:], linesty, linewidth=linewidth)
        
        if(E_fermi != None):
            axis.plot(k,np.ones(np.shape(k)) * E_fermi)
    axis.tick_params(labelsize=fontsize)
    axis.set_xticks(ticks)
    axis.set_xticklabels(label)
    axis.set_ylabel(get_unit(root, 'energy'), fontsize=int(1.5*fontsize))
    axis.set_ylim(ylim)
    axis.set_xlim([np.min(k), np.max(k)])
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


def plot_mag_cut(folder, figsize=(8,6), axis=None, fontsize=16, y_shift=0.0, color="k", dot_color="r"):
    if axis is None:
        _, axis = plt.subplots(figsize=figsize)

    cut = y_cut_idx(folder)
    x   = np.load(folder + "pos_x.npy")[cut]
    y   = np.load(folder + "pos_y.npy")[cut]
    theta = np.load(folder + "m_theta.npy")[cut]
    phi = np.load(folder + "m_phi.npy")[cut]

    m_y = np.sin(theta) * np.sin(phi)
    m_z = np.cos(theta)

    axis.quiver(-y,  x+y_shift, m_y, m_z, color=color, pivot="mid", scale=12)
    axis.plot(-y,    x+y_shift, dot_color+".")
    axis.set_aspect('equal', 'datalim')
    axis.tick_params(labelsize=fontsize)
    axis.set_yticks([])
    axis.set_xticks([])

def plot_mag_cut_v(folder, figsize=(8,6), axis=None, fontsize=16, ylim=None):
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

    axis.quiver(x, y, m_z, m_y, pivot="mid", scale=1.5, width=0.06)
    #axis.plot(x, y,"r.", markersize=3)
    axis.set_aspect('equal', 'datalim')
    axis.tick_params(labelsize=fontsize)
    axis.set_xticks([])
    axis.set_yticks([])

    if(not(ylim is None)):
        axis.set_ylim(ylim)


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
    n_loc = int(n_band /2)
    loc_dos = pdos[:n_loc,l_cut:u_cut] + pdos[n_loc:,l_cut:u_cut]

    integr = np.zeros(loc_dos.shape[0])
    for i in range(loc_dos.shape[0]):
        integr[i] = inte_dos(loc_dos[i,:], Edos)
    return integr

def round_mag(x, decimals=0):
    res = np.zeros(x.shape)
    mag    = int(np.floor(-np.log10(np.max(x))))

    return np.round(x, decimals=mag+decimals)

def plot_loc_dos(folder, E_window, figsize=(8,8), axis=None, fig=None, n_lvls=20,
        lvls=None, cmap=plt.cm.viridis, fontsize=16, lw=1):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)

    loc_dos = calc_local_dos(folder, E_window)
    x = np.load(folder + "pos_x.npy")
    y = np.load(folder + "pos_y.npy")
    if(lvls is None):
        #lvls = np.linspace(np.min(loc_dos), np.max(loc_dos), n_lvls)
        lvls = np.linspace(0.0, np.max(loc_dos), n_lvls)

    lvls = np.unique(round_mag(lvls, decimals=2))

    triang = tri.Triangulation(x, y)
    cbar = axis.tricontourf(triang, loc_dos,levels=lvls, cmap=cmap)
    axis.set_aspect('equal')

    cut = y_cut_idx(folder)
    axis.plot(x[cut], y[cut], 'C1', linewidth=lw)
    #axis.plot(x,y, 'r.', markersize=1.5)

    # axis.set_xlabel(r"$x/a_0$", fontsize=fontsize)
    # axis.set_ylabel(r"$y/a_0$", fontsize=fontsize)

    for tick in axis.xaxis.get_major_ticks():
                tick.label.set_fontsize(int(0.8*fontsize))
    for tick in axis.yaxis.get_major_ticks():
                tick.label.set_fontsize(int(0.8*fontsize))


    cb = fig.colorbar(cbar, ax=axis, ticks=lvls[::3])
    cb.ax.tick_params(labelsize=int(0.7*fontsize))

    return axis.get_xlim(), axis.get_ylim()

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

def plot_hall(folder, figsize=(8,8), axis=None, xlim=None, ylim=None, fontsize=16, linesty='', linewidth=1, label=r"$\sigma_{xy}$", color=None):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)

    E = np.load(folder + "hall_cond_E.npy")
    try:
        c = np.load(folder + "hall_cond.npy")
    except:
        cf = sorted(glob(folder + "hall_cond_iter*.npy"))[-1]
        c = np.load(cf)

    axis.plot(E,c, linesty, linewidth=linewidth, label=label, color=color)
    axis.set_ylabel(r"$\frac{\sigma_{xy}}{(e^2/h)}$", fontsize=int(1.5*fontsize))
    axis.set_xlabel(r"eV", fontsize=fontsize)
    axis.set_title(r"Hall conductance", fontsize=fontsize)

    if(not(xlim is None)):
        axis.set_xlim(xlim)
    else:
        axis.set_xlim([np.min(E), np.max(E)])
        
    if(not(ylim is None)):
        axis.set_ylim(ylim)
        
    axis.tick_params(labelsize=fontsize)
    #axis.grid(True)

def plot_hally(folder, figsize=(8,8), axis=None, xlim=None, ylim=None, fontsize=16, linesty='', linewidth=1, label="$\sigma_{xy}$ / #layers", color=None):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)

    height = lambda fol: int(np.unique(np.load(fol + "pos_z.npy")).shape[0])

    E = np.load(folder + "hall_cond_E.npy")
    try:
        c = np.load(folder + "hall_cond.npy")
    except:
        cf = sorted(glob(folder + "hall_cond_iter*.npy"))[-1]
        c = np.load(cf)
    
    c /= height(folder)

    axis.plot(E,c, linesty, linewidth=linewidth, label=label, color=color)
    axis.set_ylabel(r"$\frac{\sigma_{xy}}{(e^2/h)}$ / #layers", fontsize=int(1.5*fontsize))
    axis.set_xlabel(r"eV", fontsize=fontsize)
    axis.set_title(r"Hall conductivity", fontsize=fontsize)

    if(not(xlim is None)):
        axis.set_xlim(xlim)
    else:
        axis.set_xlim([np.min(E), np.max(E)])

    if(not(ylim is None)):
        axis.set_ylim(ylim)
        
    axis.tick_params(labelsize=fontsize)


def plot_orbmag(folder, figsize=(8,8), axis=None, xlim=None, ylim=None, labels=["L", "IC", "M"],
        fontsize=16, linesty=['','', ''], linewidth=1, which="M", color=None, prefactor=(137.03599**2)):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)

    E = np.load(folder + "orbmag_E.npy")
    if(xlim is None):
        sel = np.arange(E.shape[0])
    else:
        sel = np.argwhere(np.logical_and(E>=xlim[0], E<=xlim[1]))

    if(which.upper() == "L" or which.upper()=="ALL"):
        try:
            L = np.load(folder + "orbmag_L.npy")
        except:
            L = sorted(glob(folder + "orbmag_L_iter*.npy"))[-1]
            L = np.load(L)
        
        L *= prefactor
        axis.plot(E[sel],L[sel], linesty[1], linewidth=linewidth, label=labels[0], color=color)


    if(which.upper() == "IC" or which.upper()=="ALL"):
        try:
            IC = np.load(folder + "orbmag_IC.npy")
        except:
            IC = sorted(glob(folder + "orbmag_IC_iter*.npy"))[-1]
            IC = np.load(IC)
        IC *= prefactor
        axis.plot(E[sel],IC[sel], linesty[2], linewidth=linewidth, label=labels[1], color=color)

    
    if(which.upper() == "M" or which.upper()=="ALL"):
        try:
            om = np.load(folder + "orbmag.npy")
        except:
            om = sorted(glob(folder + "orbmag_iter*.npy"))[-1]
            om = np.load(om)
        om *= prefactor
        axis.plot(E[sel],om[sel], linesty[0], linewidth=linewidth, label=labels[2], color=color)

        
    axis.set_ylabel(r"M/$\mu_b$", fontsize=fontsize)
    axis.set_xlabel(r"eV", fontsize=fontsize)
    axis.set_title(r"Orbital magnetization", fontsize=fontsize)
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

def layer_winding(m):
    g_x = np.gradient(m, 1.0, axis=0)
    g_y = np.gradient(m, 1.0, axis=1)
    cro = np.cross(g_x, g_y)

    summe = 0.0
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            summe += np.abs(np.inner(m[i,j,:], cro[i,j,:]))
    summe *= 1.0/(4.0*np.pi)
    return summe

def winding_layer_wise(fol):
    n_atm = np.load(fol + "pos_x.npy").shape[0]

    pos = np.zeros((n_atm, 3))
    pos[:,0] = np.load(fol + "pos_x.npy")
    pos[:,1] = np.load(fol + "pos_y.npy")
    pos[:,2] = np.load(fol + "pos_z.npy")

    dim = [np.unique(pos[:,0]).shape[0],
           np.unique(pos[:,1]).shape[0],
           np.unique(pos[:,2]).shape[0]]
    
    phi   = np.load(fol + "m_phi.npy").reshape(dim,   order="F")
    theta = np.load(fol + "m_theta.npy").reshape(dim, order="F")

    m = np.zeros((dim[0], dim[1], dim[2] ,3))
    m[:,:,:,0] = np.sin(theta) * np.cos(phi)
    m[:,:,:,1] = np.sin(theta) * np.sin(phi)
    m[:,:,:,2] = np.cos(theta)

    summe = 0
    for i in range(dim[2]):
        summe += layer_winding(m[:,:,i,:])
    
    return summe

def plot_ACA(folder, figsize=(8,8), axis=None, xlim=None, ylim=None, label="ACA",
        fontsize=16, linesty=[''], linewidth=1, color=None):
    if axis is None:
        fig, axis = plt.subplots(figsize=figsize)

    E = np.load(folder + "orbmag_E.npy")
    if(xlim is None):
        sel = np.arange(E.shape[0])
    else:
        sel = np.argwhere(np.logical_and(E>=xlim[0], E<=xlim[1]))

    om = np.load(folder + "orbmag_ACA.npy")
    axis.plot(E[sel],om[sel], linesty[0], linewidth=linewidth, label=label, color=color)
    axis.set_xlabel(r"eV", fontsize=fontsize)
    axis.set_ylabel(r"$\mu_b$", fontsize=fontsize)

def plot_dbl_ldos(folder, e_range, figsize=(10,8), ax1=None, ax2=None, fig=None,
                    cmap=plt.cm.viridis, fontsize=16, lw=1):
    if((ax1 is None) or (ax2 is None) or (fig is None)):
        fig, (ax1, ax2) = plt.subplots(1,2, figsize=figsize)#,
                                       #gridspec_kw = {'width_ratios':[3, 1]})
    _, y_lim = plot_loc_dos(folder, e_range, cmap=cmap, axis=ax1, fig=fig, fontsize=fontsize, lw=lw)
    plot_mag_cut_v(folder, axis=ax2, fontsize=fontsize, ylim=y_lim)

def plot_DOS(folder, figsize=(8,8), axis=None, linewidth=2, fontsize=16, linesty=''):
    if(axis is None):
        fig, axis = plt.subplots(1,1, figsize=figsize)

    E   = np.load(folder + "DOS_E.npy")
    DOS = np.load(folder + "DOS.npy")

    axis.plot(E, DOS, linesty, linewidth=linewidth)
    axis.set_xlabel("Energy in eV", fontsize=fontsize)
    axis.set_ylabel("DOS", fontsize=fontsize)
    
    if(np.max(DOS) > 0.01):
        axis.set_ylim([-0.01, 1.1*np.max(DOS)])
    else:
        axis.set_ylim([-0.01, 0.01])

    axis.set_xlim([np.min(E), np.max(E)])
    axis.tick_params(labelsize=int(0.9*fontsize))


def plot_S(folder, figsize=(8,8), axis=None, linewidth=2, fontsize=16, linesty=''):
    if(axis is None):
        fig, axis = plt.subplots(1,1, figsize=figsize)
    
    E = np.load(folder + "spinmag_E.npy")
    S = np.load(folder + "spinmag.npy")

    axis.plot(E,S, linesty, linewidth=linewidth)
    axis.set_xlabel("Energy in eV", fontsize=fontsize)
    axis.set_ylabel("S", fontsize=fontsize)

    axis.set_xlim([np.min(E), np.max(E)])
    
    upper = np.max([0.01, 1.1*np.max(S)])
    lower = np.min([-0.01,1.1*np.min(S)])
    axis.set_ylim([lower, upper])
    axis.tick_params(labelsize=int(0.9*fontsize))


def show_image(fol, axis=None, figsize=(8,8), fontsize=16, file="render.png"):
    img=mpimg.imread(fol + file)
    if(axis is None):
        _, axis = plt.subplots(1,1, figsize=figsize)
    axis.imshow(img)
    
    axis.set_xticks([])
    axis.set_yticks([])

    axis.set_title(fol, fontsize=fontsize)


def mynote(axis, x, y, text, fontsize=16):
    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    axis.text(x, y, text, ha="center", va="center", size=fontsize, bbox=bbox_props)

def mynote_rel(axis, text, x_rel=0.05 , y_rel=0.92, fontsize=16):
    x_min, x_max = axis.get_xlim()
    y_min, y_max = axis.get_ylim()

    dx = x_max - x_min
    dy = y_max - y_min

    x = x_min + x_rel * dx
    y = y_min + y_rel * dy

    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    axis.text(x, y, text, ha="center", va="center", size=fontsize, bbox=bbox_props)


def plot_DOS_xz(folder, y_val, E_range, axis=None, fig=None,  fontsize=16, nlvls=10, color_range=None):
    E = np.load(folder + "DOS_E.npy")
    PDOS = np.load(folder + "DOS_partial.npy")
    sel = np.where(np.logical_and(E >= E_range[0], E <= E_range[1]))
    PDOS = PDOS[:,sel][:,0,:]
    PDOS = np.sum(PDOS, axis=1)

    N = int(PDOS.shape[0]/6)
    LDOS = np.zeros(N)
    for n,i in enumerate(range(0, 3*N, 3)):
        LDOS[n] = np.sum(PDOS[i:i+3]) + np.sum(PDOS[i+3*N:i+3*N+3])

    x     = np.load(folder + "pos_x.npy")
    y     = np.load(folder + "pos_y.npy")
    z     = np.load(folder + "pos_z.npy")
    theta = np.load(folder + "m_theta.npy")
    phi   = np.load(folder + "m_phi.npy")


    sel = np.abs(y-y_val) < 1e-6

    soll_shape = (30,14)

    LDOS  = LDOS[sel].reshape(soll_shape)
    x     = x[sel].reshape(soll_shape)
    y     = y[sel].reshape(soll_shape)
    z     = z[sel].reshape(soll_shape)
    theta = theta[sel].reshape(soll_shape)
    phi   = phi[sel].reshape(soll_shape)
    
    m_x  = np.sin(theta) * np.cos(phi)
    m_z  = np.cos(theta)
    if(color_range is None):
        lvls = np.linspace(np.min(LDOS), np.max(LDOS), nlvls)
    else:
        lvls = np.linspace(color_range[0], color_range[1], nlvls)

    if(axis is None):
        _, axis = plt.subplots(1,1)
    
    cbar = axis.contourf(x,z, LDOS, levels=lvls)
    axis.set_xlabel("X", fontsize=fontsize)
    axis.set_ylabel("Z", fontsize=fontsize)
    axis.set_title("Cut at y = {}".format(y_val), fontsize=fontsize)
    
    cb = fig.colorbar(cbar, ax=axis, ticks=lvls[::3])
    cb.ax.tick_params(labelsize=fontsize)

    axis.quiver(x, z, m_x, m_z, units="inches")

    axis.set_xlim=([-7,7])
    axis.set_ylim=([-0.5, 29.5])
    axis.tick_params(labelsize=fontsize)

def plot_DOS_yz(folder, x_val, E_range, axis=None, fig=None, fontsize=16, nlvls=10, color_range=None):
    E = np.load(folder + "DOS_E.npy")
    PDOS = np.load(folder + "DOS_partial.npy")
    sel = np.where(np.logical_and(E >= E_range[0], E <= E_range[1]))
    PDOS = PDOS[:,sel][:,0,:]
    PDOS = np.sum(PDOS, axis=1)

    N = int(PDOS.shape[0]/6)
    LDOS = np.zeros(N)
    for n,i in enumerate(range(0, 3*N, 3)):
        LDOS[n] = np.sum(PDOS[i:i+3]) + np.sum(PDOS[i+3*N:i+3*N+3])

    x     = np.load(folder + "pos_x.npy")
    y     = np.load(folder + "pos_y.npy")
    z     = np.load(folder + "pos_z.npy")
    theta = np.load(folder + "m_theta.npy")
    phi   = np.load(folder + "m_phi.npy")


    sel = np.abs(x-x_val) < 1e-6

    soll_shape = (30,14)

    LDOS  = LDOS[sel].reshape(soll_shape)
    x     = x[sel].reshape(soll_shape)
    y     = y[sel].reshape(soll_shape)
    z     = z[sel].reshape(soll_shape)
    theta = theta[sel].reshape(soll_shape)
    phi   = phi[sel].reshape(soll_shape)

    m_y  = np.sin(theta) * np.sin(phi)
    m_z  = np.cos(theta)
    if(color_range is None):
        lvls = np.linspace(np.min(LDOS), np.max(LDOS), nlvls)
    else:
        lvls = np.linspace(color_range[0], color_range[1], nlvls)

    lvls = np.round(lvls, decimals=3)

    if(axis is None):
        _, axis = plt.subplots(1,1)
    
    cbar = axis.contourf(y,z, LDOS, levels=lvls)
    axis.set_xlabel("Y", fontsize=fontsize)
    axis.set_ylabel("Z", fontsize=fontsize)
    axis.set_title("Cut at x = {}".format(x_val), fontsize=fontsize)

    cb = fig.colorbar(cbar, ax=axis, ticks=lvls[::3])
    cb.ax.tick_params(labelsize=fontsize)

    axis.quiver(y,z,m_y, m_z, units="inches")
    axis.set_xlim([-7,7])
    axis.set_ylim([-0.5, 29.5])
    axis.tick_params(labelsize=fontsize)


def plot_DOS_z_vert(folder, E_range, axis=None, fontsize=16, xlim=None, linewidth=1, linesty=''):
    E = np.load(folder + "DOS_E.npy")
    PDOS = np.load(folder + "DOS_partial.npy")
    sel = np.where(np.logical_and(E >= E_range[0], E <= E_range[1]))
    PDOS = PDOS[:,sel][:,0,:]
    PDOS = np.sum(PDOS, axis=1)

    N = int(PDOS.shape[0]/6)
    LDOS = np.zeros(N)
    for n,i in enumerate(range(0, 3*N, 3)):
        LDOS[n] = np.sum(PDOS[i:i+3]) + np.sum(PDOS[i+3*N:i+3*N+3])

    z     = np.load(folder + "pos_z.npy")

    soll_shape = (14,14, 30)

    LDOS  = LDOS.reshape(soll_shape,  order="F")
    z     = z.reshape(soll_shape,     order="F")

    LDOS = np.sum(LDOS, axis=(0,1))
    axis.plot(LDOS, z[0,0,:], linesty, linewidth=linewidth)

    axis.set_xlabel("LDOS", fontsize=fontsize)
    axis.set_ylabel("Z", fontsize=fontsize)

    if(not(xlim is None)):
        axis.set_xlim(xlim)

    axis.tick_params(labelsize=fontsize)