import numpy as np

# sphinx_gallery_thumbnail_number = 2
import pyvista as pv
from pyvista import examples

# Properties of SES-X air cushion
l_1 = 12  # 0.001  #  [m] length of the rectangular part of the air cushion
l_2 = 6  # 0.001  #   [m] length of the triangular part of the air cushion
b = 3.4  # [m] beam of the air cushion

# Define dimensions of the cushion area
x_s = l_1/2  # [m]
x_f = -l_1/2  # [m]
y_s = b/2  # [m]
y_p = -b/2  # [m]
x_b = -l_2 - l_1/2  # [m]


# Make data
dx = 0.1  # [m]
dy = 0.1  # [m]
x = np.arange(x_b, x_s, dx)
y = np.arange(y_p, y_s, dy)
x, y = np.meshgrid(x, y)


# Surface elevation properties
heading = 90  # [deg] wave heading relative to positive x-direction
zeta_a = 0.15  # [m] wave amplitude
wavelength = 3.83  # [m] wavelength og regular incident wave
g = 9.81  # [m/s^2] acceleration of gravity
k = 2*np.pi / wavelength  # [m^-1] wave number

beta = np.deg2rad(heading)  # [rad] wave heading in radians relative to positive x-direction

z = zeta_a * np.real(np.exp(-1j * k * x * np.cos(beta) - 1j * k * y * np.sin(beta)))

# Create and plot structured grid
n = len(z[:, 0])
m = len(z[0, :])

for i in range(n):
    for j in range(m):
        x_temp = x[i, j]
        y_temp = y[i, j]

        if x_temp < x_f and y_temp > y_s / (x_f - x_b) * (x_temp - x_b):
            z[i, j] = -99
        elif x_temp < x_f and y_temp < y_p / (x_f - x_b) * (x_temp - x_b):
            z[i, j] = -99

print(len(z[:, 0]))
print(len(z[0, :]))

'''
grid = pv.StructuredGrid(x, y, z)
grid.plot()
'''

x = x.flatten()
y = y.flatten()
z = z.flatten()
index_remove = np.where(z != -99)

x = x[index_remove]
y = y[index_remove]
z = z[index_remove]

points = np.c_[x, y, z]
cloud = pv.PolyData(points)
cloud.plot(point_size=2)
surf = cloud.delaunay_2d()
surf.plot(show_edges=True)
