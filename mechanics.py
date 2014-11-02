import numpy as np

class Body(object):
	"""Represents a physical body"""
	def __init__(self, name, mass):
		self.mass = mass
		self.name = name
		self.connections = set()

	def __repr__(self):
		return "<Body(%r, mass=%.3f) [conns=%d]>" % (
			self.name,
			self.mass,
			len(self.connections)
		)

class Conn(object):
	"""Represents a physical connection"""
	def __init__(self, k, lam=0):
		self.k = k
		self.lam = lam

		self.a = None
		self.b = None

	def between(self, a, b):
		self.a = a
		self.b = b

		a.connections.add(self)
		b.connections.add(self)

		return self

	def __repr__(self):
		return "<Connection(k=%r, lam=%r) between %r and %r>" % (
			self.k,
			self.lam,
			self.a.name,
			self.b.name
		)

class System(object):
	def __init__(self, containing):
		self.bodies = []
		self.connections = set()

		self.add_connected(containing)

	def add_connected(self, body):
		""" Add everything attached to $body to this system """
		if body in self.bodies:
			return
		else:
			self.bodies.append(body)
			for connection in body.connections:
				if connection.a == body:
					other = connection.b
				else:
					other = connection.a

				self.connections.add(connection)

				self.add_connected(other)

		self.bodies = [body for body in self.bodies if body.mass != float('inf')]

	@property
	def DOF(self):
		""" degrees of freedom of the system """
		return len(self.bodies)

	@property
	def mass_matrix(self):
		return np.array([body.mass for body in self.bodies])

	@property
	def k_and_lam_matrices(self):
		N = self.DOF
		k_matrix = np.zeros((N, N))
		lam_matrix = np.zeros((N, N))

		for connection in self.connections:
			# convert bodies into a number, that corresponds to their index
			# in the mass matrix
			a_is_gnd = b_is_gnd = False
			try:
				a_i = self.bodies.index(connection.a)
			except ValueError:
				a_is_gnd = True

			try:
				b_i = self.bodies.index(connection.b)
			except ValueError:
				b_is_gnd = True


			if a_is_gnd and b_is_gnd:
				raise ValueError("You're a moron")
			elif a_is_gnd:
				k_matrix[b_i, b_i] += connection.k
				lam_matrix[b_i, b_i] += connection.lam
			elif b_is_gnd:
				k_matrix[a_i, a_i] += connection.k
				lam_matrix[a_i, a_i] += connection.lam
			else:
				# force is proportional to displacement/velocity of this end
				k_matrix[a_i, a_i] += connection.k
				lam_matrix[a_i, a_i] += connection.lam
				k_matrix[b_i, b_i] += connection.k
				lam_matrix[b_i, b_i] += connection.lam

				# ... minus that of the other end
				k_matrix[a_i, b_i] -= connection.k
				lam_matrix[a_i, b_i] -= connection.lam
				k_matrix[b_i, a_i] -= connection.k
				lam_matrix[b_i, a_i] -= connection.lam

		return k_matrix, lam_matrix

	def idx(self, body):
		return self.bodies.index(body)

	def __getitem__(self, name):
		item = next(body for body in self.bodies if body.name == name)
		if item is not None:
			return item
		else:
			raise KeyError(name)

	def _unpack(self, nparray):
		""" convert a 2d array into a dictionary of 1d arrays """
		return {
			body: nparray[:,i]
			for i, body in enumerate(self.bodies)
		}



def simulate(system, forces, dt):
	"""
	Simulate system under $forces[body] sampled every $dt seconds

	"""
	steps = len(forces.values()[0])

	# convert force dict into force vector
	f = np.zeros((steps, system.DOF))
	for body, force in forces.items():
		f[:, system.idx(body)] = force

	# system properties
	m      = system.mass_matrix
	k, lam = system.k_and_lam_matrices

	# transient state
	y         = np.zeros(system.DOF)
	y_dot     = np.zeros(system.DOF)
	y_dot_dot = np.zeros(system.DOF)
	y_prev = y

	# saved state
	y_values = []
	yd_values = []
	ydd_values = []

	# do the main integration loop
	for i in range(0, steps):
		# total acceleration
		y_dot_dot = (f[i] - k.dot(y) - lam.dot(y_dot))/m

		# verlet integrate
		y_prev, y = y, (2*y - y_prev + dt*dt*y_dot_dot)

		# estimate velocity
		y_dot = (y - y_prev)/dt

		y_values.append(y)
		yd_values.append(y_dot)
		ydd_values.append(y_dot_dot)


	y_values = np.array(y_values)
	yd_values = np.array(yd_values)
	ydd_values = np.array(ydd_values)

	return system._unpack(y_values), system._unpack(yd_values), system._unpack(ydd_values)

def frequency_response(system, omegas, shape):
	m      = np.diag(system.mass_matrix)
	k, lam = system.k_and_lam_matrices

	# convert force dict into force vector
	f = np.zeros(system.DOF)
	for body, force in shape.items():
		f[system.idx(body)] = force

	results = []
	for w in omegas:
		matrix = (m * (1.j * w) ** 2 + lam * (1.j * w) + k)
		result = np.dot(np.linalg.inv(matrix), f)
		results.append(result)


	return system._unpack(
		np.array(results)
	)
