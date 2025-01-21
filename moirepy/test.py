import numpy as np
import matplotlib.pyplot as plt

# Define the rhombus parameters
lv1 = np.array([1, -0.3])
lv2 = np.array([-0.5, np.sqrt(3)/2])
mln1 = 1
mln2 = 1

class Rhombus:
    def __init__(self):
        self.mlv1 = lv1
        self.mlv2 = lv2                     
        self.mln1 = mln1
        self.mln2 = mln2

    def inside_rhombus(self, points):
        p, q = points.T

        a1, b1 = self.mln1 * self.mlv1
        a2, b2 = self.mln2 * self.mlv2

        m1 = b1 / a1
        c1 = b1 - m1*a1
        c2 = b2 + b1 - m1*(a1+a2)
        if c1 > c2: c1, c2 = c2, c1

        m2 = b2 / a2
        d1 = b2 - m2*a2
        d2 = b1 + b2 - m2*(a1+a2)
        if d1 > d2: d1, d2 = d2, d1

        return (c1 <= q - m1*p) & (q - m1*p < c2) & (d1 <= q - m2*p) & (q - m2*p < d2)

# Generate random points
np.random.seed(0)
all_points = (mln1 + mln2) * (2*np.random.rand(10000, 2) - 1)

# Create an instance of the Rhombus class
rhombus = Rhombus()

# Filter the points using the inside_rhombus function
inside_mask = rhombus.inside_rhombus(all_points)
inside_points = all_points[inside_mask]

# Plot both sets of points
plt.scatter(all_points[:, 0], all_points[:, 1], c='red', s=50, alpha=0.5)
plt.scatter(inside_points[:, 0], inside_points[:, 1], c='green', s=10, alpha=1)


# plot: lv1*mln1, lv2*mln2, origin and lv1*mln1+lv2*mln2
plt.scatter(
	[0, lv1[0]*mln1, lv2[0]*mln2, lv1[0]*mln1 + lv2[0]*mln2],
	[0, lv1[1]*mln1, lv2[1]*mln2, lv1[1]*mln1 + lv2[1]*mln2],
	c='k',
	s=50,
	alpha=1
)


plt.axis('equal')
plt.show()
