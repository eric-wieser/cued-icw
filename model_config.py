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

	import json

	with open('optimized.json') as f:
		absorber_params = json.load(f)

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
