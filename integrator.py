import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.signal

import model_config
import inputs

import mechanics

system = model_config.system

dt = 1/2000.0
t = np.arange(0, 30, dt)
f = scipy.signal.chirp(t, 1, 30, 20)
freqs = scipy.interp(t, [0, 30], [1, 20])

y, yd, ydd = mechanics.simulate(
	system,
	forces={
		model_config.floors[0]: f
	},
	dt=dt
)

# do displacement plots
_, (disp_ax, ab_disp_ax, force_ax) = plt.subplots(3, 1, sharex=True)

for floor in model_config.floors:
	disp_ax.plot(freqs, y[floor], label=floor.name, linewidth=0.5)

disp_ax.set_title('output displacement')
disp_ax.grid()
disp_ax.legend()

for ab in model_config.absorbers:
	ab_disp_ax.plot(freqs, y[ab], label=ab.name, linewidth=0.5)

ab_disp_ax.set_title('output displacement')
ab_disp_ax.grid()
ab_disp_ax.legend()

force_ax.plot(freqs, f, linewidth=0.5)
force_ax.set_title('input force')
force_ax.grid()


# Do fft plots
def fourier(signal):
	phasors = np.fft.fft(signal)
	phasors = phasors[:len(phasors)/2]

	freqs = np.fft.fftfreq(signal.size, d=dt)
	freqs = freqs[:len(freqs)/2]

	return freqs, np.abs(phasors)

_, (disp_fft_ax, ab_disp_fft_ax, force_fft_ax) = plt.subplots(3, 1, sharex=True)

for floor in model_config.floors:
	disp_fft_ax.plot(
		*fourier(y[floor]),
		label=floor.name
	, linewidth=0.5)
disp_fft_ax.set_xlim(0, 30)
disp_fft_ax.set_title('output fft')
disp_fft_ax.grid()
disp_fft_ax.legend()

for ab in model_config.absorbers:
	i = system.idx(ab)
	ab_disp_fft_ax.plot(
		*fourier(y[ab]),
		label=ab.name,
		linewidth=0.5)
ab_disp_fft_ax.set_xlim(0, 30)
ab_disp_fft_ax.set_title('input fft')
ab_disp_fft_ax.grid()
ab_disp_fft_ax.legend()

force_fft_ax.plot(*fourier(f), linewidth=0.5)
force_fft_ax.set_xlim(0, 30)
force_fft_ax.set_title('input fft')
force_fft_ax.grid()

plt.show()
