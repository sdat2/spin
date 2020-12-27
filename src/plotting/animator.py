#!/usr/bin/env python3
"""animator.py by sdat2
This program is designed to take the pickled output
and plot them to produce an effective animation. """
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter
import src.plotting.my_plotting_style as mps
import src.time_wrapper as twr
import src.data_loading.io_pickles as ip
import src.plotting.plotters as plot
import src.simulation.sclass as scl
Writer = animation.writers["ffmpeg"]
writer = Writer(fps=80, metadata=dict(artist="sdat2"), bitrate=1800)


class AniMP4:
    """A class to animate things"""

    @twr.timeit
    def __init__(
        self, OUTPUT_DIREC, name="forgot_name_TC_Animator", move_with=True, **kwargs
    ):
        """load the pickle files into memory,
        then graphing objects to self"""
        self.co, self.sy = ip.read_in(OUTPUT_DIREC, name)
        # set up figure and animation
        self.move_with = move_with
        self.name = name
        #  Replot Energies etc.
        with plt.style.context(("fivethirtyeight")):  # plot the normal things
            plt.clf()
            plot_energy(self.co, self.sy)
            plot_angular_momentum(self.co, self.sy)
            plot_multi_AM(self.co, self.sy)
        with plt.style.context(("dark_background")):
            self.fig = plt.figure()  # the main animation figure
            self.ax = self.fig.add_subplot(
                111, aspect="equal", autoscale_on=False, xlim=(-20, 20), ylim=(-20, 20)
            )  # ax has all these useful attributes
            #  self.ax.grid() # Add gridlines
            self.no_particles = len(self.sy.coordinate_grid[0, :, 0])
            if self.no_particles < 10:  # Deal with multi body animate test case
                (self.line,) = self.ax.plot(
                    [], [], "wo", markersize=1, label="Test Masses"
                )
                (self.galactic_centre,) = self.ax.plot(
                    [0], [0], "wo", markersize=1, label="Galaxy A"
                )
                (self.impactor,) = self.ax.plot(
                    [0], [0], "wo", markersize=1, label="Galaxy B"
                )
            else:
                (self.line,) = self.ax.plot(
                    [], [], "wo", markersize=0.5, label="Test Masses"
                )
                (self.galactic_centre,) = self.ax.plot(
                    [0], [0], color="red", marker="+", label="Galaxy A"
                )
                (self.impactor,) = self.ax.plot(
                    [0], [0], color="green", marker="+", label="Galaxy B"
                )
            if self.no_particles > 10:
                if self.co.halo == True:
                    if self.move_with == True:
                        print("I want to plot a Halo.")
            print("May have plotted Halo if needed")

            if move_with:  # In the non inertial coordinate case
                plt.xlabel(r"$X^{\prime}$")
                plt.ylabel(r"$Y^{\prime}$")
            else:
                plt.xlabel(r"$X$")
                plt.ylabel(r"$Y$")

            self.add_details()  # Comment out line if this is TMI

            # Add some text outputs in suitable areas of the figure
            self.time_text = self.ax.text(1.01, 0.95, "", transform=self.ax.transAxes)
            self.KE_text = self.ax.text(1.01, 0.88, "", transform=self.ax.transAxes)
            self.GPE_text = self.ax.text(1.01, 0.81, "", transform=self.ax.transAxes)
            self.energy_text = self.ax.text(1.01, 0.74, "", transform=self.ax.transAxes)
            plt.tight_layout()  # This supposedly makes stops the label from falling off.

        print("Timer is " + str(len(self.sy.short_timer)) + " Long")
        self.dt = (
            self.sy.timer[1] - self.sy.timer[0]
        )  # read the time step directly from the timer file

    def add_details(self):
        """Add the details of which algorithm timestep etc, to make the
        animations more informative on there own"""

        if self.co.algorithm == "vv":
            algo = "Verlocity Verlot"
        if self.co.algorithm == "rk4o":
            algo = "Runge Kutta Forth Order"
        if self.co.algorithm == "herm":
            algo = "Hermite Fourth Order"

        self.algorithm_title = self.ax.text(
            1.01, 0.65, "Algorithm:", transform=self.ax.transAxes
        )
        self.algorithm_text = self.ax.text(
            1.01, 0.58, algo, transform=self.ax.transAxes
        )
        self.timestep_text = self.ax.text(
            1.01, 0.51, "dt =" + str(self.co.tstep), transform=self.ax.transAxes
        )
        self.length_softening_distance = self.ax.text(
            1.01,
            0.44,
            r"$\epsilon$ = " + str(self.co.epsilon),
            transform=self.ax.transAxes,
        )

    def ani_init(self):
        """initialize animation with lots of blanks where data will go"""
        self.line.set_data([], [])
        self.galactic_centre.set_data([], [])
        self.impactor.set_data([], [])
        self.time_text.set_text("")
        self.KE_text.set_text("")
        self.GPE_text.set_text("")
        self.energy_text.set_text("")
        return (
            self.line,
            self.galactic_centre,
            self.impactor,
            self.time_text,
            self.KE_text,
            self.GPE_text,
            self.energy_text,
        )
        # One might hope that there was a better way to return objects

    def animate(self, i):
        """perform an animation step in the comoving or inertial frame"""
        # save the cooridnates of the first particle as local vaiables
        self.time_text.set_text("Time = %.1f" % self.sy.short_timer[i])
        self.KE_text.set_text("KE = %.6f" % self.sy.Total_Energies["Kinetic"][i])
        self.GPE_text.set_text(
            "GPE = %.6f" % self.sy.Total_Energies["Gravitational"][i]
        )
        self.energy_text.set_text("E = %.6f" % self.sy.Total_Energies["Total"][i])
        plt.tight_layout()  # This supposedly makes stops the label from falling off.
        a = self.sy.coordinate_grid[0, 0, i]
        b = self.sy.coordinate_grid[1, 0, i]
        if self.move_with == False:  # inertial frame
            self.line.set_data(
                self.sy.coordinate_grid[0, :, i], self.sy.coordinate_grid[1, :, i]
            )
            self.galactic_centre.set_data(a, b)
            self.impactor.set_data(
                self.sy.coordinate_grid[0, -1, i], self.sy.coordinate_grid[1, -1, i]
            )
        if self.move_with == True:  # comoving frame
            self.line.set_data(
                self.sy.coordinate_grid[0, :, i] - a,
                self.sy.coordinate_grid[1, :, i] - b,
            )
            self.galactic_centre.set_data(0.0, 0.0)
            self.impactor.set_data(
                self.sy.coordinate_grid[0, -1, i] - a,
                self.sy.coordinate_grid[1, -1, i] - b,
            )
        return (
            self.line,
            self.galactic_centre,
            self.impactor,
            self.time_text,
            self.KE_text,
            self.GPE_text,
            self.energy_text,
        )

    @twr.timeit
    def animate_starter(self, **kwargs):
        """ Get the animation itself started """
        interval = 5  # this number works fine, but is rather arbirtary, presumably in milliseconds
        print("The timer length is " + str(len(self.sy.short_timer)))
        print("Shape of coordinate_grid is " + str(np.shape(self.sy.coordinate_grid)))
        print("The animator interval was " + str(interval) + " in unknown units")
        # I don't currently understand why the galaxy chooses
        # to slow down mid way through.
        # Perhaps I should look at the FuncAnimation
        #  dictionary and work out what has gone wrong.
        with plt.style.context(("dark_background")):
            ani = animation.FuncAnimation(
                self.fig,
                self.animate,
                frames=len(self.sy.short_timer),
                interval=interval,
                blit=True,
                init_func=self.ani_init,
            )
            ani.save(
                str(self.co.out)
                + "/"
                + str(self.name)
                + "move_with_"
                + str(self.move_with)
                + ".mp4",
                writer=writer,
            )
            plt.clf()  # always make sure you close the lid
