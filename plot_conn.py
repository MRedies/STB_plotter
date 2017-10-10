import numpy as np
from mayavi.mlab import *
from colorsys import hsv_to_rgb
import sys
import os.path
import glob


def plot_conn(file):
    x         = np.load(file + "pos_x.npy")
    y         = np.load(file + "pos_y.npy")
    z         = np.load(file + "pos_z.npy")
    neigh     = np.load(file + "neigh.npy") - 1
    conn_type = np.load(file + "conn_type.npy")
    print neigh

    l = points3d(x, y, z, resolution=50, scale_factor=0.3, color=(0.2, 0.4, 0.5))


    for atom in range(neigh.shape[0]):
        for n in range(neigh.shape[1]):
            if(neigh[atom,n] >= 0 and conn_type[atom,n] == 1): 
                px = [x[atom], x[n]]
                print len(px)
                py = [y[atom], y[n]]
                pz = [z[atom], z[n]]

                plot3d(px,py,pz)#, tube_radius=2)
                break

plot_conn("/home/matthias/STB/output/dbg/")
