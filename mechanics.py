import numpy as np

class Body(object):
	"""Represents a physical body"""
	def __init__(self, mass, name):
		self.mass = mass
		self.name = name
		self.connections = set()

	def __repr__(self):
		return "<Body %r with mass %.3f>"%(self.name, self.mass)

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
		return "<Connection between %r and %r>"%(self.a.name, self.b.name)
	
Ground = Body(float('inf'), 'ground')


f1 = Body(1.83, name='Floor 1')
f2 = Body(1.83, name='Floor 2')
f3 = Body(1.83, name='Floor 3')
a1 = Body(0.183, name='Absorber 1')
a2 = Body(0.183, name='Absorber 2')

Conn(k=4200,lam=1.98).between(Ground, f1)
Conn(k=4200,lam=1.98).between(f1, f2)
Conn(k=4200,lam=1.98).between(f2, f3)
Conn(k=585.22,lam=1.98).between(f3, a1)
Conn(k=1220.9,lam=1.98).between(f3, a2)

def simulate(body):
	bodies = set()
	connections = set()

	def find_all(body):
		if body in bodies:
			return
		else:
			bodies.add(body)
			for connection in body.connections:
				if connection.a == body:
					other = connection.b
				else:
					other = connection.a

				connections.add(connection)

				find_all(other)

	find_all(body)

	print bodies
	print connections

	bodies = list(body for body in bodies if body.mass != float('inf'))
	N = len(bodies)

	mass_matrix = np.array([body.mass for body in bodies])

	k_matrix = np.zeros((N, N))
	lam_matrix = np.zeros((N, N))

	for connection in connections:

		# if either end is connected to ground
		if connection.a.mass == float('inf'):
			b_i = bodies.index(connection.b)
			k_matrix[b_i, b_i] += connection.k
			lam_matrix[b_i, b_i] += connection.lam

		elif connection.b.mass == float('inf'):
			a_i = bodies.index(connection.a)
			k_matrix[a_i, a_i] += connection.k
			lam_matrix[a_i, a_i] += connection.lam

		else:
			a_i = bodies.index(connection.a)
			b_i = bodies.index(connection.b)

			# (k) effect is proportional to displacement of this end
			k_matrix[a_i, a_i] += connection.k
			lam_matrix[a_i, a_i] += connection.lam
			k_matrix[b_i, b_i] += connection.k
			lam_matrix[b_i, b_i] += connection.lam


			# ... minus displacement of other
			k_matrix[a_i, b_i] -= connection.k
			lam_matrix[a_i, b_i] -= connection.lam
			k_matrix[b_i, a_i] -= connection.k
			lam_matrix[b_i, a_i] -= connection.lam


	return mass_matrix, k_matrix, lam_matrix, bodies
