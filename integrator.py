import matplotlib.pyplot as plt
import numpy as np
import math

import model_config
import inputs

N = 3
t = 0
y_dot = np.zeros(N)
y = np.zeros(N)
dt = 0.0005
max_t = 30

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


lam = make_mat({
	(None, 0): 1.98,
	(0,    1): 1.98,
	(1,    2): 1.98
})

# lam = make_mat({})
m = model_config.mass



def fourier(signal):
	phasors = np.fft.fft(signal)
	phasors = phasors[:len(phasors)/2]

	freqs = np.fft.fftfreq(signal.size, d=dt)
	freqs = freqs[:len(freqs)/2]

	return freqs, np.abs(phasors)


# f = f_sin_maker(13.73912302)
# f = inputs.on_floor(inputs.pulse, floor=0, n=N)
f = inputs.on_floor(inputs.from_data(1), floor=0, n=N)

y_prev = y
y_values = []
f_values = []
t_values = []

while t < max_t:
	y_dot_dot = (f(t) - k.dot(y) - lam.dot(y_dot))/m
	temp = y
	y = 2*y - y_prev + dt*dt*y_dot_dot
	y_prev = temp
	y_dot = (y - y_prev)/dt

	y_values.append(y)
	f_values.append(f(t))
	t_values.append(t)

	t = t+dt


y_values = np.array(y_values)
f_values = np.array(f_values)

fig, (disp_ax, freq_ax, force_ax) = plt.subplots(3, 1)

for i in range(N):
	disp_ax.plot(t_values, y_values[:,i], label = "floor %d" % i)
disp_ax.set_title('displacement')
disp_ax.grid()
disp_ax.legend()

for i in range(N):
	freq_ax.plot(
		*fourier(y_values[:,i]),
		label="floor %d" % i
	)
freq_ax.set_xlim(0, 30)
freq_ax.set_title('frequency')
freq_ax.grid()
freq_ax.legend()


for i in range(N):
	force_ax.plot(t_values, f_values[:,i], label = "floor %d" % i)
force_ax.set_title('force')
force_ax.grid()
force_ax.legend()

plt.show()
