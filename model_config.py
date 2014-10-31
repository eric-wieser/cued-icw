# physical parameters
mass = 1.83  # of one floor
L = 0.2      # height of wall
E = 210E9    # Young's Modulus
b = 0.08     # width of wall plate
d = 0.001    # thickness of wall plate

# calculated properties
I = b*d*d*d/12  # second moment of area
k = (24*E*I)/(L*L*L)  # static stiffness for each floor


from mechanics import Body, Conn, System
import math

def make_building():
	ground = Body('ground', mass=float('inf'))
	f1 = Body('Floor 1',    mass=mass)
	f2 = Body('Floor 2',    mass=mass)
	f3 = Body('Floor 3',    mass=mass)
	floors = [f1, f2, f3]

	Conn(k=k, lam=1.98).between(ground, f1)
	Conn(k=k, lam=1.98).between(f1, f2)
	Conn(k=k, lam=1.98).between(f2, f3)

	return ground, floors


def make_absorber(freq, attached_to, lam=1.98, mass=0.1*mass):
	a = Body('%.2fHz on %s' % (freq, attached_to.name), mass=mass)
	omega = 2 * math.pi * freq
	c = Conn(k=omega**2 * mass, lam=lam).between(a, attached_to)
	print c.k
	return a
