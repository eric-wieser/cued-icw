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

	@property
	def DOF(self):
		""" degrees of freedom of the system (sort of) """
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
			a_i = self.bodies.index(connection.a)
			b_i = self.bodies.index(connection.b)

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
