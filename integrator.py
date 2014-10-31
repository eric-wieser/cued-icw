import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.signal

import model_config
import inputs
from mechanics import Conn, Body, System, simulate, frequency_response

# configure input signals
dt = 1/100.0 # dt reduced for efficiency in looping process
t = np.arange(0, 30, dt)
f = scipy.signal.chirp(t, 1, 30, 20)
freqs = scipy.interp(t, [0, 30], [1, 20])
omegas = freqs * np.pi * 2

print omegas

# perform iterative process to distribute mass over many dampers

no_iterations = 10 # no dampers we wish to add
res_freq_values = []
total_damper_mass = model_config.mass


def do_plots():

	# do displacement plots
	_, (disp_ax, ab_disp_ax, force_ax) = plt.subplots(3, 1, sharex=True)

	for floor in floors:
		disp_ax.plot(freqs, y[floor], label=floor.name, linewidth=0.5)

	disp_ax.set_title('output displacement')
	disp_ax.grid()
	disp_ax.legend()

	for ab in absorbers:
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

	_, (disp_fft_ax, disp_fft2_ax, ab_disp_fft_ax, force_fft_ax) = plt.subplots(4, 1, sharex=True)

	for floor in floors:
		disp_fft_ax.plot(
			*fourier(y[floor]),
			label=floor.name
		, linewidth=0.5)
	disp_fft_ax.set_xlim(0, 30)
	disp_fft_ax.set_title('output fft')
	disp_fft_ax.grid()
	disp_fft_ax.legend()

	resp = frequency_response(system, omegas=omegas, shape={ floors[0]: 1 })


	for floor in floors:
		print resp[floor].shape
		disp_fft2_ax.plot(
			freqs,
			abs(resp[floor]),
			label=floor.name
		, linewidth=0.5)
	disp_fft2_ax.set_xlim(0, 30)
	disp_fft2_ax.set_title('output fft')
	disp_fft2_ax.grid()
	disp_fft2_ax.legend()

	for ab in absorbers:
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


first_time = True

while len(res_freq_values) <= no_iterations:
	print "current absorbers: %r" % res_freq_values

	# take the base building, and add absorber
	ground, floors = model_config.make_building()
	absorbers = []

	for fr in res_freq_values:
		a = model_config.make_absorber(freq=fr, attached_to=floors[2], mass=total_damper_mass/len(res_freq_values))
		absorbers.append(a)
	system = System(containing=ground)

	# simulate
	y, yd, ydd = simulate(
		system,
		forces={
			floors[0]: f
		},
		dt=dt
	)

	if first_time:
		do_plots()
		first_time = False

	# get peak frequency of max amplitude and add to array
	top_floor = floors[-1]
	max_i = y[top_floor].argmax()
	max_f = freqs[max_i]
	res_freq_values.append(max_f)

do_plots()

plt.show()
