import numpy as np
import scipy as sc
import scipy.linalg
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)

# physical parameters
N = 3  # number of degrees of freedom
mass = 1.83  # of one floor
L = 0.2      # height of wall
E = 210E9    # Young's Modulus
b = 0.08     # width of wall plate
d = 0.001    # thickness of wall plate

# calculated properties
I = b*d*d*d/12  # second moment of area
k = (24*E*I)/(L*L*L)  # static stiffness for each floor

# matrix representation
# To include vibration absorbers, you will need to modify these
M = np.diag([mass, mass, mass])
K = k*np.array([
	[ 2, -1,  0],
	[-1,  2, -1],
	[ 0, -1,  1]
])


eig_values, eig_vects = sc.linalg.eig(K, M)

freqs = np.sqrt(eig_values)
modeshapes = eig_vects

#  Print natural frequencies and mode vectors in command window
print "Hertz:", freqs/(2*np.pi)
print "Shapes:", modeshapes

# harmonic solution for unit force at floor 1
omegas = np.r_[0:130:0.01]
displ = np.array([
	np.dot(
		np.linalg.inv(K - (w ** 2) * M),
		[1, 0, 0]
	)
	for w in omegas
])

# guess discontinuities and mark them as infinite
displ[displ > 0.01] = float('inf')
displ[displ < -0.01] = -float('inf')


def linear_plot():
	fig, ax = plt.subplots(1, 1)
	ax.plot(omegas, displ[:,0], color='blue',  label='floor 1')
	ax.plot(omegas, displ[:,1], color='green', label='floor 2')
	ax.plot(omegas, displ[:,2], color='red',   label='floor 3')
	ax.set_xlim(omegas[0], omegas[-1])
	ax.set_xlabel(R'\omega')
	ax.set_ylim(-1e-3, 1e-3)
	ax.set_ylabel('displacement / m')
	ax.grid()
	ax.legend()


def log_plot():
	fig, ax = plt.subplots(1, 1)
	ax.plot(omegas, abs(displ[:,0]), color='blue',  label='floor 1')
	ax.plot(omegas, abs(displ[:,1]), color='green', label='floor 2')
	ax.plot(omegas, abs(displ[:,2]), color='red',   label='floor 3')
	ax.set_xlim(omegas[0], omegas[-1])
	ax.set_xlabel(R'\omega')
	ax.set_yscale('log')
	ax.set_ylabel('displacement / m')
	ax.grid()
	ax.legend()


def plot_modeshapes():
	# add in a zero point for the ground floor
	centers = np.vstack([
		[0, 0, 0],
		modeshapes
	])

	# find coordinates of left and right wall vertices
	left = centers - 0.25
	right = centers + 0.25
	heights = np.r_[0:N+1]

	fig, axes = plt.subplots(1, N)
	for i, ax in enumerate(axes):
		ax.axis([-2, 2, 0, 3.5])
		ax.grid()
		ax.set_title('Mode {}'.format(i))

		# draw walls
		ax.plot(left[:,i], heights, color='black')
		ax.plot(right[:,i], heights, color='black')

		# draw floors
		for x1, x2, y in zip(left[:,i], right[:,i], heights):
			ax.plot([x1, x2], [y, y], color='black')


linear_plot()
log_plot()
plot_modeshapes()

plt.show()