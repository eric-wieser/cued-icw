import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.signal

import model_config
import inputs
from mechanics import Conn, Body, System, simulate, frequency_response

# configure input signals
dt = 1/200.0 # dt reduced for efficiency in looping process
t = np.arange(0, 90, dt)
freqs = scipy.interp(t, [0, t.max()], [0, 20])
omegas = freqs * np.pi * 2

# perform iterative process to distribute mass over many dampers

no_iterations = 25 # no dampers we wish to add
absorber_params = []
total_damper_mass = model_config.mass # can be varied, chosen arbitrarily

do_optimize = False

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
max_amps = []

def optimize():
	while len(absorber_params) <= no_iterations:
		print "Step %d/%d" % (len(absorber_params), no_iterations)
		# take the base building
		ground, floors = model_config.make_building()
		absorbers = []

		for fr, attached_to in absorber_params:
			a = model_config.make_absorber(
				freq=fr,
				attached_to=floors[attached_to],
				mass=total_damper_mass/len(absorber_params),
				lam=model_config.lam/len(absorber_params)
			)
			absorbers.append(a)
		system = System(containing=ground)

		freq_resp = frequency_response(system, omegas=omegas, shape={ floors[0]: 1 })

		fig, freq_plot = plt.subplots()
		fig.figurePatch.set_alpha(0)
		for floor, fcolor in zip(floors, floor_colors):
			freq_plot.plot(freqs, np.abs(freq_resp[floor]), color=fcolor, linewidth=0.5)
		freq_plot.set(
			xlabel="Frequency / Hz",
			ylabel="Amplitude / m",
			ylim=[0, 0.005],
			title="{} absorber{}".format(
				len(absorber_params),
				's' if len(absorber_params) != 1 else ' '
			)
		)
		freq_plot.grid()
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
		
		max_amps.append(max_amp)

		absorber_params.append((max_f, max_floor_i))

if do_optimize:
	optimize()

	for frequency, floor in absorber_params:
		print "floor: %d \t freq: %.2f" % (floor, frequency)

	# plot max amp vs no dampers
	fig, axis = plt.subplots()
	fig.figurePatch.set_alpha(0)
	axis.bar(np.arange(len(max_amps)), max_amps)
	axis.grid()
	axis.set(xlabel="number of absorbers", ylabel="maximum harmonic response (all floors)")
	fig.savefig("graphs/amp-vs-no-abs.png")
	plt.show()

else:
	ground, floors = model_config.make_optimized_building()
	system = System(ground)

# simulate
s_ground, s_floors = model_config.make_building()
s_system = System(s_ground)

input_data = [
	('step', inputs.step(t, at=10)),
	('sweep', scipy.signal.chirp(t, 0, t.max(), 20)),
	('earthquake', inputs.from_data(t, 3))
]

shared_fig, axes = plt.subplots(len(input_data), sharex=True)

for shared_ax, (name, f) in zip(axes, input_data):
	y_before, _, _ = simulate(
		s_system,
		forces={
			s_floors[0]: f
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

	fig, (before_dplot, after_dplot, f_plot) = plt.subplots(3, 1, sharex=True)
	fig.figurePatch.set_alpha(0)

	max_i_disp = np.max(y_before.values())
	min_i_disp = np.min(y_before.values())
	i_disp_range = max_i_disp - min_i_disp
	ylim = [min_i_disp - 0.1*i_disp_range, max_i_disp + 0.1*i_disp_range]

	# plot floors before
	for floor, fcolor in reversed(zip(s_floors, floor_colors)):
		before_dplot.plot(t, y_before[floor], label=floor.name, color=fcolor, linewidth=0.5)
	before_dplot.set(
		ylim=ylim,
		title="without absorbers",
		ylabel="displacement / m"
	)
	before_dplot.grid()
	before_dplot.legend()

	# plot floors after
	for floor, fcolor in reversed(zip(floors, floor_colors)):
		after_dplot.plot(t, y_after[floor], label=floor.name, color=fcolor, linewidth=0.5)
	after_dplot.set(
		ylim=ylim,
		title="with absorbers",
		ylabel="displacement / m"
	)
	after_dplot.grid()

	# plot forces
	f_range = f.max() - f.min()
	f_plot.plot(t, f, linewidth=0.5, color='k')
	f_plot.set(
		ylabel='input force / N',
		xlabel='time / s',
		ylim=[f.min() - 0.1*f_range, f.max() + 0.1*f_range]
	)
	f_plot.grid()

	fig.savefig('graphs/disp-{}.png'.format(name))

	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
	ax2.plot(*fourier(f), color="k", linewidth=0.5)
	ax2.set(
		ylabel='amplitude / N',
		xlabel='frequency / Hz',
		xlim=(0,20)
	)

	resp = frequency_response(s_system, omegas, { s_floors[0]: 1 })

	for floor, fcolor in reversed(zip(s_floors, floor_colors)):
		ax1.plot(freqs, np.abs(resp[floor]), label=floor.name, color=fcolor, linewidth=0.5)
	ax1.set(
		ylabel='amplitude / N',
		xlabel='frequency / Hz'
	)

	fig.savefig('graphs/fft-{}.png'.format(name))

	for floor, fcolor in reversed(zip(s_floors, floor_colors)):
		shared_ax.plot(t, y_before[floor], label=floor.name, color=fcolor, linewidth=0.5)

	shared_ax.set(title=name)


shared_fig.savefig('graphs/disp-all.png'.format(name))