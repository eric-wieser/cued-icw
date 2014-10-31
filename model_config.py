# physical parameters
mass = 1.83  # of one floor
L = 0.2      # height of wall
E = 210E9    # Young's Modulus
b = 0.08     # width of wall plate
d = 0.001    # thickness of wall plate

# calculated properties
I = b*d*d*d/12  # second moment of area
k = (24*E*I)/(L*L*L)  # static stiffness for each floor
