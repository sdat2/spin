import numpy as np
import matplotlib.pyplot as plt
import src.plotting.my_plotting_style as mps
import src.simulation.halos as hal
import src.simulation.sclass as scl


def acc_fun_plotter():
    """ A function designed to plot the gravitational field around a central particle"""
    co = scl.Controls(EPS=0.01)
    x_values = np.linspace(0, 25, num=2000)
    # x_lim = [-0.1 10]
    x_lim = [-0.1, 25]
    a_values = []
    # radii = [0.1, 0.2, 0.3, 0.5]
    radii = [6, 7, 8, 0.1]
    for radius_index in range(len(radii)):
        a_values.append([])
        test_particle = [
            hal.Particle(
                np.array([0, 0, 0]),
                np.array([0, 0, 0]),
                1.0,
                True,
                0,
                halo=True,
                radius=radii[radius_index],
            )
        ]

        for x in x_values:
            pos = np.array([x, 0, 0])
            acc_vec = -hal.acc_fun(-1, pos, test_particle, co)
            a_values[radius_index].append(acc_vec[0])

        plt.plot(
            x_values,
            a_values[radius_index],
            LineWidth=0.5,
            label="Halo with radius " + str(radii[radius_index]),
        )
        
    y_lim = [-0.005, 1.4 * (np.max(a_values[0]))]

    plt.ylim(y_lim)
    plt.gca().invert_xaxis()  # this command seems to do absolutely nothing :/
    plt.ylabel(r"Attractive Acceleration")
    plt.xlabel(r"Separation Distance")
    plt.xlim(x_lim)
    plt.legend()
    plt.tight_layout()
    # plt.savefig('Gen_Output/Test_Halo_Forces.pdf')
    plt.savefig("Gen_Output/Dark_vs_Light.pdf")
    plt.clf()


if __name__ == "__main__":
    with plt.style.context(("fivethirtyeight")):
        acc_fun_plotter()
