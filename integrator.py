import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.signal

import model_config
import inputs
from mechanics import Conn, Body, System, simulate, frequency_response

# configure input signals
dt = 1/200.0 # dt reduced for efficiency in looping process
t = np.arange(0, 30, dt)
f = scipy.signal.chirp(t, 1, 30, 20)
freqs = scipy.interp(t, [0, 30], [1, 20])
omegas = freqs * np.pi * 2

print omegas

# perform iterative process to distribute mass over many dampers

no_iterations = 50 # no dampers we wish to add
absorber_params = []
total_damper_mass = model_config.mass # can be varied, chosen arbitrarily

floor_colors = ["r", "g", "b"]

# do displacement plots
def do_disp_plots():
	_, (disp_ax, ab_disp_ax, force_ax) = plt.subplots(3, 1, sharex=True)

	for floor, fcolor in zip(floors, floor_colors):
		disp_ax.plot(freqs, y[floor], label=floor.name, color=fcolor, linewidth=0.5)

	disp_ax.set_title('output displacement')
	disp_ax.grid()
	disp_ax.legend()

	for ab in absorbers:
		floor_id = floors.index(ab.floor)
		ab_disp_ax.plot(freqs, y[ab], label=ab.name, color=floor_colors[floor_id], linewidth=0.5)

	ab_disp_ax.set_title('output displacement')
	ab_disp_ax.grid()
	ab_disp_ax.legend(loc=2,prop={'size':6})

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

def do_dfft_plots():
	_, (disp_fft_ax, ab_disp_fft_ax, force_fft_ax) = plt.subplots(3, 1, sharex=True)

	for floor, fcolor in zip(floors, floor_colors):
		disp_fft_ax.plot(
			*fourier(y[floor]),
			label=floor.name,
			color=fcolor
		, linewidth=0.5)
	disp_fft_ax.set_xlim(0, 30)
	disp_fft_ax.set_title('output fft')
	disp_fft_ax.grid()
	disp_fft_ax.legend()

	for ab in absorbers:
		i = system.idx(ab)
		floor_id = floors.index(ab.floor)
		ab_disp_fft_ax.plot(
			*fourier(y[ab]),
			color=floor_colors[floor_id], 
			label=ab.name,
			linewidth=0.5
		)
	ab_disp_fft_ax.set_xlim(0, 30)
	ab_disp_fft_ax.set_title('input fft')
	ab_disp_fft_ax.grid()

	force_fft_ax.plot(*fourier(f), linewidth=0.5)
	force_fft_ax.set_xlim(0, 30)
	force_fft_ax.set_title('input fft')
	force_fft_ax.grid()

def do_plots():
	do_disp_plots()
	do_dfft_plots()
	do_freq_plot()

def do_freq_plot():
	_, ax = plt.subplots()

	for floor, fcolor in zip(floors, floor_colors):
		ax.plot(
			freqs,
			abs(freq_resp[floor]),
			label=floor.name,
			color=fcolor,
			linewidth=0.5
		)
	ax.set_xlim(0, 30)
	ax.set_title('output fft')
	ax.grid()
	ax.legend()

first_time = True

while len(absorber_params) <= no_iterations:
	print "Step %d/%d" % (len(absorber_params), no_iterations)
	# take the base building
	ground, floors = model_config.make_building()
	absorbers = []

	for fr, attached_to in absorber_params:
		a = model_config.make_absorber(
			freq=fr,
			attached_to=floors[attached_to],
			mass=total_damper_mass/len(absorber_params)
		)
		absorbers.append(a)
	system = System(containing=ground)

	freq_resp = frequency_response(system, omegas=omegas, shape={ floors[0]: 1 })

	fig, freq_plot = plt.subplots()
	fig.figurePatch.set_alpha(0)
	for floor, fcolor in zip(floors, floor_colors):
		freq_plot.plot(freqs, np.abs(freq_resp[floor]), color=fcolor, linewidth=0.5)
	freq_plot.set_xlabel("Frequency / Hz")
	freq_plot.set_ylabel("Amplitude / m")
	freq_plot.set_ylim(0, 0.020)
	fig.savefig('graphs/absorber-{:02d}.png'.format(len(absorber_params)))
	plt.close(fig)

	max_by_floor = [
		(
			i,
			freqs[np.abs(freq_resp[floor]).argmax()],
			np.abs(freq_resp[floor]).max()
		)
		for i, floor in enumerate(floors)
	]
	max_floor_i, max_f, max_amp = max(
		max_by_floor,
		key=lambda (floor, freq, amp): amp
	)
	
	absorber_params.append((max_f, max_floor_i))

for frequency, floor in absorber_params:
	print "floor: %d \t freq: %.2f" % (floor, frequency)

# simulate
s_ground, s_floors = model_config.make_building()
s_system = System(s_ground)

y_before, _, _ = simulate(
	s_system,
	forces={
		floors[0]: f
	},
	dt=dt
)
y_after, _, _ = simulate(
	system,
	forces={
		floors[0]: f
	},
	dt=dt
)

_, (before_dplot, after_dplot, f_plot) = plt.subplots(3, 1, sharex=True)

# plot floors before
for floor, fcolor in zip(s_floors, floor_colors):
	before_dplot.plot(freqs, y_before[floor], label=floor.name, color=fcolor, linewidth=0.5)
before_dplot.set_title('output displacement')
before_dplot.grid()
before_dplot.legend()

# plot floors after
for floor, fcolor in zip(floors, floor_colors):
	after_dplot.plot(freqs, y_after[floor], label=floor.name, color=fcolor, linewidth=0.5)

after_dplot.set_title('output displacement')
after_dplot.grid()
after_dplot.legend()

# plot forces
f_plot.plot(freqs, f, linewidth=0.5)
f_plot.set_title('input force')
f_plot.grid()


do_plots()

plt.show()
