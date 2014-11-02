import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import scipy.io

import model_config
import inputs
from mechanics import System, simulate

ground, floors = model_config.make_building()
system = System(ground)

# load in the real data
data = [
	sc.io.loadmat('f-floor0.mat', squeeze_me=True),
	sc.io.loadmat('f-floor1.mat', squeeze_me=True),
	sc.io.loadmat('f-floor2.mat', squeeze_me=True)
]
f = np.mean([d['data_ch1'] for d in data], axis=0)
t = data[0]['t']
floor_d = [d['data_ch2'] for d in data]


_, _, ydd = simulate(system, { floors[0]: f}, np.diff(t)[0])

fig, (ax3, ax1, ax2) = plt.subplots(3, 1, sharex=True)

ax1.set(title="measured response", ylabel="acceleration")
ax2.set(title="simulated response", ylabel=u"acceleration")
ax3.set(title="input", ylabel="force ")

legends = []

for i, floor in enumerate(floors):
	ax1.plot(t, floor_d[i], label=floor.name, color=floor.color, linewidth=0.5)
	l, = ax2.plot(t, ydd[floor], label=floor.name, color=floor.color, linewidth=0.5)
	legends.append((l, floor.name))


ax1.grid()
ax2.grid()

ax3.plot(t, f)
ax3.grid()

lines, labels = zip(*legends)
fig.legend(lines, labels, 'lower center', ncol=3)


plt.show()