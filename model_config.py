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

ground = Body('ground', mass=float('inf'))
f1 = Body('Floor 1',    mass=mass)
f2 = Body('Floor 2',    mass=mass)
f3 = Body('Floor 3',    mass=mass)
floors = [f1, f2, f3]

a1 = Body('Absorber 1', mass=mass/10)
a2 = Body('Absorber 2', mass=mass/10)
absorbers = [a1, a2]

Conn(k=k,      lam=1.98).between(ground, f1)
Conn(k=k,      lam=1.98).between(f1, f2)
Conn(k=k,      lam=1.98).between(f2, f3)
Conn(k=585.22, lam=1.98).between(f3, a1)
Conn(k=1220.9, lam=1.98).between(f3, a2)

system = System(containing=ground)
