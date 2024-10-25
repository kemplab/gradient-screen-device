import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # appropriate import to draw 3d polygons
from matplotlib import style
import matplotlib

font = {'family' : 'normal',
        'weight' : 'regular',
        'size'   : 12}

matplotlib.rc('font', **font)

plt.figure('SPLTV',figsize=(10,5))
ax = plt.subplot(121, projection='3d')


def plot_tern_plane(ax, vertices):
    a = vertices[0]
    b = vertices[1]
    c = vertices[2]

    x1 = np.array([a, 0, 0])
    y1 = np.array([0, b, 0])
    z1 = np.array([0, 0, c])
    #ax.scatter(x1, y1, z1)

    verts = [list(zip(x1, y1, z1))]
    srf = Poly3DCollection(verts, alpha = 0.40, facecolor = 'crimson')
    ax.add_collection3d(srf)



plot_tern_plane(ax, (1, 1, 1))
plot_tern_plane(ax, (2, 2, 2))
plot_tern_plane(ax, (4, 4, 4))
plot_tern_plane(ax, (8, 8, 8))
#plot_tern_plane(ax, (16, 16, 16))

#x = np.array([0, .01, .1, 10])
x = np.array([0, 1, 2, 4, 8])
y = x
z = x


xgrid, ygrid, zgrid = np.meshgrid(x, y, z)

#ax = plt.subplot(111, projection = "3d")
ax.scatter(xgrid, ygrid, zgrid, c = "blue", label = "Well Plate")

# axis specific stuff
ax.set_xlabel('Drug 1')
ax.set_ylabel('Drug 2')
ax.set_zlabel('Drug 3')

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# Now set color to white (or whatever is "invisible")
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')


#ax.legend()
ax.grid(False)

plt.show()