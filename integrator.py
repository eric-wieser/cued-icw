import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.signal

import model_config
import inputs
from mechanics import System, simulate, frequency_response

# Sampling rate and range
dt = 1/200.0
t = np.arange(0, 90, dt)

# frequency sweeping
freqs = scipy.interp(t, [0, t.max()], [0, 20])
omegas = freqs * np.pi * 2

do_optimize = True

def fourier(signal):
	phasors = np.fft.fft(signal)
	phasors = phasors[:len(phasors)/2]

	freqs = np.fft.fftfreq(signal.size, d=dt)
	freqs = freqs[:len(freqs)/2]

	return freqs, np.abs(phasors)

def optimize():
	"""perform iterative process to distribute mass over many dampers"""

	max_amps = []
	no_iterations = 25 # no dampers we wish to add
	total_damper_mass = model_config.mass # can be varied, chosen arbitrarily
	absorber_params = []

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
		for floor in floors:
			freq_plot.plot(freqs, np.abs(freq_resp[floor]), color=floor.color, linewidth=0.5)
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


	# plot max amp vs no dampers
	fig, axis = plt.subplots()
	fig.figurePatch.set_alpha(0)
	axis.bar(np.arange(len(max_amps)), max_amps)
	axis.grid()
	axis.set(xlabel="number of absorbers", ylabel="maximum harmonic response (all floors)")
	fig.savefig("graphs/amp-vs-no-abs.png")
	plt.show()

	return absorber_params

if do_optimize:
	absorber_params = optimize()

	import json
	with open('optimized.json', 'w') as f:
		json.dump(absorber_params, f)

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

	# fiddle with scales
	def padded_lim(d, p=0.1):
		max_d = np.max(d)
		min_d = np.min(d)
		range_d = max_d - min_d

		return [min_d - range_d*p, max_d + range_d*p]

	ylim = padded_lim(y_before.values())

	# plot floors before addition of absorbers
	for floor in reversed(s_floors):
		before_dplot.plot(t, y_before[floor], label=floor.name, color=floor.color, linewidth=0.5)
	before_dplot.set(
		ylim=ylim,
		title="without absorbers",
		ylabel="displacement / m"
	)
	before_dplot.grid()
	before_dplot.legend()

	# plot floors after addition of absorbers
	for floor in reversed(floors):
		after_dplot.plot(t, y_after[floor], label=floor.name, color=floor.color, linewidth=0.5)
	after_dplot.set(
		ylim=ylim,
		title="with absorbers",
		ylabel="displacement / m"
	)
	after_dplot.grid()

	# plot forces
	f_plot.plot(t, f, linewidth=0.5, color='k')
	f_plot.set(
		ylabel='input force / N',
		xlabel='time / s',
		ylim=padded_lim(f)
	)
	f_plot.grid()

	fig.savefig('graphs/disp-{}.png'.format(name))


	fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)

	# plot frequency spectrum of earthquake
	ax2.plot(*fourier(f), color="k", linewidth=0.5)
	ax2.set(
		ylabel='amplitude / N',
		xlabel='frequency / Hz',
		xlim=(0,20)
	)

	# plot frequency response of initial system
	resp = frequency_response(s_system, omegas, { s_floors[0]: 1 })
	for floor in reversed(s_floors):
		ax1.plot(freqs, np.abs(resp[floor]),
			label=floor.name, color=floor.color, linewidth=0.5
		)
	ax1.set(
		ylabel='amplitude / N',
		xlabel='frequency / Hz'
	)

	fig.savefig('graphs/fft-{}.png'.format(name))


	# add this displacement plot to the shared graph
	for floor in reversed(s_floors):
		shared_ax.plot(t, y_before[floor], label=floor.name, color=floor.color, linewidth=0.5)
	shared_ax.set(title=name)


shared_fig.savefig('graphs/disp-all.png'.format(name))
