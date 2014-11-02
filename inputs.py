import numpy as np
import scipy.io

def sine(t, freq_hz):
	return np.sin(2 * np.pi * freq_hz * t)


def step(t, at=1):
	return np.where(t > at, 1, 0)

def from_data(t, file_id):
	fnames = {
		1: 'dataset1-small.mat',
		2: 'dataset2-moderate.mat',
		3: 'dataset3-large.mat',
	}
	fname = fnames[file_id]

	f_force = scipy.io.loadmat(fname, squeeze_me=True)['f0']
	f_t = np.arange(len(f_force)) * 1.0/500  # 500 Hz
	f_t += 1 # wait 1 second before earthquake begins

	if f_t[-1] > t[-1]:
		import warnings
		warnings.warn(
			"Earthquake continues beyond sampling range (until {}s)".format(f_t[-1]),
			UserWarning
		)

	return np.interp(t, f_t, f_force, left=0, right=0)

