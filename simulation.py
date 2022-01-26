import pandas as pd
import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


plt.close('all')

n_particles = 686

times = []

particles = {new_list: [] for new_list in range(n_particles)}

n = 0

df = pd.read_csv("galaxies.data",
                 delim_whitespace=True,
                 usecols=[1,2,3],
                 names=["x", "y", "z"])
for i in range(len(df)):
    x = df.loc[i,"x"]
    y = df.loc[i, "y"]
    z = df.loc[i, "z"]
    n = (n + 1) % n_particles
    particles[n].append((x,y,z))

categories = np.zeros(n_particles, dtype=np.int8)
for i in range(len(categories)):
    categories[i] = int(i % 2)
colormap = np.array(['cyan', 'red'])


xdata, ydata, zdata = [], [], []
for n in range(n_particles):
        x = particles[n][i][0]
        y = particles[n][i][1]
        z = particles[n][i][2]
        xdata.append(x)
        ydata.append(y)
        zdata.append(z)

plt.style.use("dark_background")
fig = plt.figure(figsize=(20,20))
ax = fig.add_subplot(projection='3d')
graph = ax.scatter(xdata, ydata, zdata, s=8, color=colormap[categories])

# set plot details
axis_lim = 120
ax.view_init(90,90)

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('w')
ax.yaxis.pane.set_edgecolor('w')
ax.zaxis.pane.set_edgecolor('w')

ax.set_xlim(-axis_lim, axis_lim)
ax.set_ylim(-axis_lim, axis_lim)
ax.set_zlim(-axis_lim, axis_lim)
ax.set_xlabel("x (kpc)")
ax.set_ylabel("y (kpc)")
ax.set_zlabel("z (kpc)")
plt.axis('off')
ax.set_aspect("auto")
ax.grid(False)


def animate(i):
    xdata, ydata, zdata = [], [], []
    for n in range(n_particles):
        x = particles[n][i][0]
        y = particles[n][i][1]
        z = particles[n][i][2]
        xdata.append(x)
        ydata.append(y)
        zdata.append(z)
    npx = np.array(xdata)
    npy = np.array(ydata)
    npz = np.array(zdata)

    graph._offsets3d = (npx, npy, npz)


anim = animation.FuncAnimation(fig, animate, frames=333, interval=10, blit=False)

anim.save('simulation.mp4', writer='ffmpeg', fps=24)
plt.show()