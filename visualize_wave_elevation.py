import numpy as np
import pyvista as pv


def create_unstructured_mesh(heading, frequency, t=0, zeta_a=0.2, dx=0.1, dy=0.1, l_1=12, l_2=6, b=3.4):

    # Define dimensions of the cushion area
    x_s = l_1 / 2  # [m]
    x_f = -l_1 / 2  # [m]
    y_s = b / 2  # [m]
    y_p = -b / 2  # [m]
    x_b = -l_2 - l_1 / 2  # [m]

    # Make data
    x = np.arange(x_b, x_s, dx)
    y = np.arange(y_p, y_s, dy)
    x, y = np.meshgrid(x, y)

    beta = np.deg2rad(heading)  # [rad] wave heading in radians relative to positive x-direction
    k = frequency**2/g

    z = zeta_a * np.real(np.exp(1j * k * (frequency * t - x * np.cos(beta) - y * np.sin(beta))))

    n = len(z[:, 0])  # number of rows
    m = len(z[0, :])  # number of columns

    # Remove grid points from
    for i in range(n):
        for j in range(m):
            if x[i, j] < x_f and (y[i, j] > y_s / (x_f - x_b) * (x[i, j] - x_b) or
                                  y[i, j] < y_p / (x_f - x_b) * (x[i, j] - x_b)):
                z[i, j] = -99

    x = x.flatten()
    y = y.flatten()
    z = z.flatten()
    index_remove = np.where(z != -99)

    points = np.c_[x[index_remove], y[index_remove], z[index_remove]]

    return points


def set_visualization_settings():
    '''Defines visualization settings to be used when plotting in Pyvista

    :return:
    '''
    # Set white background color and black font of labels to fit with document style.
    pv.set_plot_theme("document")
    # Center scalar bar.
    pv.global_theme.colorbar_horizontal.position_x = 0.5 - pv.global_theme.colorbar_horizontal.width / 2


if __name__ == "__main__":
    set_visualization_settings()

    # Refinement of the mesh
    dx = 0.05  # [m]
    dy = 0.05  # [m]

    heading = 0  # [deg] wave heading relative to positive x-direction
    zeta_a = 0.1  # [m] wave amplitude
    g = 9.81  # [m/s^2] acceleration of gravity
    wavelength = 10.3  # [m] wavelength og regular incident wave
    k = 2*np.pi / wavelength  # [m^-1] wave number
    frequency = np.sqrt(g * k)  # [rad/s] incident wave frequency
    t = 0  #0.10*2*np.pi/frequency  #  [s] time instance

    points = create_unstructured_mesh(heading, frequency, t, zeta_a, dx, dy)

    cloud = pv.PolyData(points)
    surf = cloud.delaunay_2d()
    plotter = pv.Plotter(shape=(1, 1))


    plotter.add_mesh(surf, scalars=surf.points[:, 2], cmap="CET_L6")

    plotter.view_xy()
    #plotter.camera.focal_point = (0.2, 0.3, 0.3)
    #plotter.camera.up = (0.0, 1.0, 0.0)
    plotter.camera.zoom(1.2)
    plotter.show()
