import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.signal
import model_config
import inputs

N = 5
t = 0
y_dot = np.zeros(N)
y = np.zeros(N)
dt = 1/2000.0
max_t = 30
no_floors = 3

def make_mat(conns):
	mat = np.zeros((N,N))

	def get_conn(i,j):
		return conns.get((i,j), conns.get((j,i), 0))

	# get non-diagonal matrix values. k[i,j] = k[j,i] = -(spring connecting i with j)
	for i in range(N):
		for j in range(N):
			if i > j:
				mat[i,j] = mat[j,i] = -get_conn(i,j)

	for i in range(N):
		mat[i,i] = sum(
			get_conn(i,k)
			for k in range(N)
		) + get_conn(i,None)

	return mat

k = model_config.k * make_mat({
	(None, 0): 1,
	(0,    1): 1,
	(1,    2): 1
})

k += make_mat({
	(2,    3): 585.22,
	(2,    4): 1220.9
})


lam = make_mat({
	(None, 0): 1.98,
	(0,    1): 1.98,
	(1,    2): 1.98,
	(2,    3): 1.98, ##
	(2,    4): 1.98 ##
})

names = [""] * N
names[:no_floors] = ('floor %d'    % i for i in range(no_floors))
names[no_floors:] = ('absorber %d' % i for i in range(N - no_floors))

# lam = make_mat({})
m = np.zeros(N)
m[:no_floors] = model_config.mass
m[no_floors:] = model_config.mass/10

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

f = inputs.on_floor(f, 0, N)
# f = inputs.on_floor(inputs.pulse, floor=0, n=N)
# f = inputs.on_floor(inputs.from_data(3), floor=0, n=N)

y_prev = y
y_values = []
f_values = []
t_values = []
a_values = []
freq_values = []

while t < max_t:

	y_dot_dot = (f(t) - k.dot(y) - lam.dot(y_dot))/m
	#y_dot_dot[3] = y_dot_dot[3] * 10; ##
	#y_dot_dot[4] = y_dot_dot[4] * 10; ##
	temp = y
	y = 2*y - y_prev + dt*dt*y_dot_dot
	y_prev = temp
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

for i in range(no_floors):
	disp_ax.plot(freq_values, y_values[:,i], label=names[i], linewidth=0.5)
disp_ax.set_title('output displacement')
disp_ax.grid()
disp_ax.legend()

for i in range(no_floors, N):
	ab_disp_ax.plot(freq_values, y_values[:,i], label=names[i], linewidth=0.5)
ab_disp_ax.set_title('output displacement')
ab_disp_ax.grid()
ab_disp_ax.legend()

for i in range(no_floors):
	force_ax.plot(freq_values, f_values[:,i], label=names[i], linewidth=0.5)
force_ax.set_title('input force')
force_ax.grid()
force_ax.legend()

# Do fft plots
_, (disp_fft_ax, ab_disp_fft_ax, force_fft_ax) = plt.subplots(3, 1, sharex=True)

for i in range(no_floors):
	disp_fft_ax.plot(
		*fourier(y_values[:,i]),
		label=names[i]
	, linewidth=0.5)
disp_fft_ax.set_xlim(0, 30)
disp_fft_ax.set_title('output fft')
disp_fft_ax.grid()
disp_fft_ax.legend()

for i in range(no_floors, N):
	ab_disp_fft_ax.plot(
		*fourier(y_values[:,i]),
		label=names[i]
	, linewidth=0.5)
ab_disp_fft_ax.set_xlim(0, 30)
ab_disp_fft_ax.set_title('input fft')
ab_disp_fft_ax.grid()
ab_disp_fft_ax.legend()

for i in range(no_floors):
	force_fft_ax.plot(
		*fourier(f_values[:,i]),
		label=names[i]
	, linewidth=0.5)
force_fft_ax.set_xlim(0, 30)
force_fft_ax.set_title('input fft')
force_fft_ax.grid()
force_fft_ax.legend()

plt.show()
