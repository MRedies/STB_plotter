import numpy as np
from mayavi.mlab import *
from colorsys import hsv_to_rgb
import sys
import os.path
import glob

def color(theta, phi):
    while(phi < 0):
        phi += 2.0*np.pi
    h = phi / (2.0 * np.pi)
    s = 1
    v = 1.0 - 0.9999999*(theta / np.pi)
    print("h = {}, s = {}, v ={}".format(h,s,v))
    return hsv_to_rgb(h,s,v)


def plot_mag(file):
    x = np.load(file + "pos_x.npy")
    y = np.load(file + "pos_y.npy")
    z = np.load(file + "pos_z.npy")

    phi   = np.load(file + 'm_phi.npy')
    theta = np.load(file + 'm_theta.npy')

    u = np.sin(theta) * np.cos(phi)
    v = np.sin(theta) * np.sin(phi)
    w = np.cos(theta)

    obj = quiver3d(x, y, z, u, v, w,
            line_width=3, colormap='hsv',
            scale_factor=3, mode='arrow',resolution=25)
    # for i in range(x.shape[0]):
        # r,g,b = color(theta[i], phi[i])
        # print("R: {}, G: {}, B: {}".format(r,g,b))
        # obj = quiver3d(x[i], y[i], z[i], u[i], v[i], w[i],
                # line_width=3, color=(r,g,b), colormap='hsv',
                # scale_factor=0.8, mode='arrow',resolution=25)
    return obj

files = sorted(glob.glob("/home/matthias/Masterthesis/figs/insulator_skyrm/insulator_band_large_lin_skyrm_apd=*/"))
for i in files:
    if(not os.path.isfile(i + "/skyrm.png")):
        figure(bgcolor=(1,1,1))
        plot_mag(i + "/")
        savefig(i + "/skyrm.png")
        #close()
