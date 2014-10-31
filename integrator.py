import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.signal

import model_config
import inputs

system = model_config.system
m      = system.mass_matrix
k, lam = system.k_and_lam_matrices
N      = system.DOF

t = 0
y_dot = np.zeros(N)
y = np.zeros(N)
dt = 1/2000.0
max_t = 30
no_floors = 3


def fourier(signal):
	phasors = np.fft.fft(signal)
	phasors = phasors[:len(phasors)/2]

	freqs = np.fft.fftfreq(signal.size, d=dt)
	freqs = freqs[:len(freqs)/2]

	return freqs, np.abs(phasors)


#f = inputs.sin_maker(1, 10)

def f(t):
	return scipy.signal.chirp(t, 1, 30, 20)

def freq(t):
	return scipy.interp(t, [0, 30], [1, 20])

f = inputs.on_floor(f, system.idx(model_config.floors[0]), N)

y_prev = y
y_values = []
f_values = []
t_values = []
a_values = []
freq_values = []

# do the main integration loop
while t < max_t:
	# total acceleration
	y_dot_dot = (f(t) - k.dot(y) - lam.dot(y_dot))/m

	# verlet integrate
	y_prev, y = y, (2*y - y_prev + dt*dt*y_dot_dot)

	# estimate velocity
	y_dot = (y - y_prev)/dt

	y_values.append(y)
	f_values.append(f(t))
	t_values.append(t)
	a_values.append(y_dot_dot)
	freq_values.append(freq(t))

	t = t+dt


y_values = np.array(y_values)
f_values = np.array(f_values)
a_values = np.array(a_values)
freq_values = np.array(freq_values)

# do displacement plots
_, (disp_ax, ab_disp_ax, force_ax) = plt.subplots(3, 1, sharex=True)

for floor in model_config.floors:
	i = system.idx(floor)
	disp_ax.plot(freq_values, y_values[:,i], label=floor.name, linewidth=0.5)

disp_ax.set_title('output displacement')
disp_ax.grid()
disp_ax.legend()

for ab in model_config.absorbers:
	i = system.idx(ab)
	ab_disp_ax.plot(freq_values, y_values[:,i], label=ab.name, linewidth=0.5)

ab_disp_ax.set_title('output displacement')
ab_disp_ax.grid()
ab_disp_ax.legend()

for i, body in enumerate(system.bodies):
	if not np.all(f_values[:,i] == 0):
		force_ax.plot(freq_values, f_values[:,i], label=body.name, linewidth=0.5)
force_ax.set_title('input force')
force_ax.grid()
force_ax.legend()

# Do fft plots
_, (disp_fft_ax, ab_disp_fft_ax, force_fft_ax) = plt.subplots(3, 1, sharex=True)

for floor in model_config.floors:
	i = system.idx(floor)
	disp_fft_ax.plot(
		*fourier(y_values[:,i]),
		label=floor.name
	, linewidth=0.5)
disp_fft_ax.set_xlim(0, 30)
disp_fft_ax.set_title('output fft')
disp_fft_ax.grid()
disp_fft_ax.legend()

for ab in model_config.absorbers:
	i = system.idx(ab)
	ab_disp_fft_ax.plot(
		*fourier(y_values[:,i]),
		label=ab.name
	, linewidth=0.5)
ab_disp_fft_ax.set_xlim(0, 30)
ab_disp_fft_ax.set_title('input fft')
ab_disp_fft_ax.grid()
ab_disp_fft_ax.legend()

for i, body in enumerate(system.bodies):
	if not np.all(f_values[:,i] == 0):
		force_fft_ax.plot(
			*fourier(f_values[:,i]),
			label=body.name
		, linewidth=0.5)
force_fft_ax.set_xlim(0, 30)
force_fft_ax.set_title('input fft')
force_fft_ax.grid()
force_fft_ax.legend()

plt.show()
