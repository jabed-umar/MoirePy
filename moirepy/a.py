import numpy as np
import matplotlib.pyplot as plt
from layers import TriangularLayer

lattice1 = TriangularLayer()
lattice1.generate_points2(7*2+1)

lattice2 = TriangularLayer()
# lattice2.perform_rotation(16.26020471)
lattice2.perform_rotation(9.43)
lattice2.generate_points2(7*2+1)

plt.plot(lattice1.points[:,0], lattice1.points[:,1], 'o')
plt.plot(lattice2.points[:,0], lattice2.points[:,1], 'o')


plt.axvline(x=0, color='k')
plt.axhline(y=0, color='k')

plt.show()