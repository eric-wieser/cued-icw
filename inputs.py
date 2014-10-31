import math
import numpy as np
import scipy.io

def sin_maker(amp, freq_hz):
	def sin(t):
		return amp*math.sin(2 * math.pi * freq_hz * t)

	return sin

def pulse(t):
	if t >= 1:
		return 1
	return 0

def from_array(a, sample_rate):
	def f(t):
		t_i = (t - 1) // 0.002
		if t_i < len(force):
			return force[t_i]

		# after sample
		return None

	return f


def on_floor(func, floor, n):
	def f(t):
		res = np.zeros(n)
		res[floor] = func(t)
		return res

	return f

def from_data(file_id):
	fnames = {
		1: 'dataset1-small.mat',
		2: 'dataset2-moderate.mat',
		3: 'dataset3-large.mat',
	}
	fname = fnames[file_id]

	all_data = scipy.io.loadmat(fname, squeeze_me=True)
	force = all_data['f0']

	def f(t):
		# before earthquake
		if t < 1:
			return 0

		t_i = (t - 1) // 0.002
		if t_i < len(force):
			return force[t_i]

		# after earthquake
		return 0

	return f

