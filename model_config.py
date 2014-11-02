# physical parameters
mass = 1.83  # of one floor
L = 0.2      # height of wall
E = 210E9    # Young's Modulus
b = 0.08     # width of wall plate
d = 0.001    # thickness of wall plate
lam = 1.98

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

	f1.color = 'r'
	f2.color = 'g'
	f3.color = 'b'

	Conn(k=k, lam=lam).between(ground, f1)
	Conn(k=k, lam=lam).between(f1, f2)
	Conn(k=k, lam=lam).between(f2, f3)

	return ground, floors


def make_absorber(freq, attached_to, lam=lam, mass=0.1*mass):
	a = Body('%.2fHz on %s' % (freq, attached_to.name), mass=mass)
	omega = 2 * math.pi * freq

	k = omega**2 * mass # + lam**2 / (4 * mass)
	c = Conn(k=k, lam=lam).between(a, attached_to)

	a.floor = attached_to

	return a

def make_optimized_building():
	total_damper_mass = mass

	absorber_params = [
		(3.39, 2),
		(2.35, 2),
		(4.48, 2),
		(2.11, 2),
		(9.74, 0),
		(2.93, 2),
		(8.37, 0),
		(3.79, 2),
		(2.68, 2),
		(4.82, 1),
		(7.78, 2),
		(7.43, 0),
		(3.14, 2),
		(2.54, 2),
		(4.00, 2),
		(2.26, 2),
		(2.04, 2),
		(10.84, 0),
		(1.98, 2),
		(3.55, 2),
		(9.18, 0),
		(2.82, 2),
		(4.17, 2),
		(4.99, 2),
		(2.48, 2),
		(5.18, 2)
	]

	ground, floors = make_building()
	absorbers = []

	for fr, attached_to in absorber_params:
		a = make_absorber(
			freq=fr,
			attached_to=floors[attached_to],
			mass=total_damper_mass/len(absorber_params),
			lam=lam/len(absorber_params)
		)
		absorbers.append(a)

	return ground, floors
