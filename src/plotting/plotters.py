#!/usr/bin/env python3
'''plotters.py by sdat2
This program is designed to take the pickled output
and plot them to produce an effective animation '''
import os # needed for making dircetories.
import numpy as np
import matplotlib.pyplot as plt
from my_plotting_style import *
from sclass import *
import os


def plot_energy(co, sy):
    ''' Plot Energies'''
    print('plotting energies with the name ' + co.name)
    print('giving it the directory with name ' + co.out)
    plt.plot(sy.short_timer, sy.Total_Energies['Kinetic'], label='KE', LineWidth=1)
    plt.plot(sy.short_timer, sy.Total_Energies['Gravitational'], label='GPE', LineWidth=1)
    plt.plot(sy.short_timer, sy.Total_Energies['Total'], label='Total Energy', LineWidth=1)
    x_lim = [np.min(sy.short_timer), np.max(sy.short_timer)]
    plt.xlim(x_lim) # Don't waste length.
    plt.plot(x_lim, [0,0], color='k', LineWidth=1)  # Put an origin line in
    # Plot the other energy types
    plt.xlabel(r'Simulation Time')
    # Units for both of these quantities are unclear
    plt.ylabel(r'Simulation Energy')
    plt.legend()  # plot the legend
    plt.tight_layout()  # should keep the figure bounds ok
    plt.savefig(str(co.out)+'/'+str(co.name)+'_Energies_Equal.pdf')
    plt.clf()  # clear the figure

def plot_angular_momentum(co, sy):
    ''' Plot Angular Momentum '''
    print('plotting angular momentum with the name ' + co.name)
    print('giving it the directory with name ' + co.out)
    plt.plot(sy.short_timer, sy.Angular_Momentum['Total'][2], label='Total AM', LineWidth=1)
    x_lim = [np.min(sy.short_timer), np.max(sy.short_timer)]
    plt.xlim(x_lim) # Don't waste length.
    plt.plot(x_lim, [0,0], color='k', LineWidth=1)  # Put an origin line in
    plt.xlabel(r'Simulation Time')  # Units for both of these quantities are unclear
    plt.ylabel(r'Angular Momentum')
    plt.legend()
    plt.tight_layout()
    plt.savefig(co.out+'/'+co.name+'_Total_AM.pdf')
    plt.clf()

def plot_multi_AM(co, sy):
    '''Plot many angular momentums'''
    print('plotting angular momentums with the name ' + co.name)
    print('giving it the directory with name ' + co.out)
    #no_particles = len(sy.coordinate_grid[0,:,0])
    for i in sy.Angular_Momentum:
        plt.plot(sy.short_timer, sy.Angular_Momentum[i][2], label='AM '+str(i), LineWidth=1)
    x_lim = [np.min(sy.short_timer), np.max(sy.short_timer)]
    plt.xlim(x_lim) # Don't waste length.
    plt.plot(x_lim, [0,0], color='k', LineWidth=1)  # Put an origin line in
    plt.xlabel(r'Simulation Time')  # Units for both of these quantities are unclear
    plt.ylabel(r'Angular Momentum')
    plt.legend()
    plt.tight_layout()
    plt.savefig(co.out+'/'+co.name+'_AMs.pdf')
    plt.clf()  # clear the figure


@timeit
def snap_shot(co, particles, move_with=True, Time=-1000):
    '''Take a snap shot at a point in the simulation.
    Possibly include some Velocity information'''
    print('Taking a snapshot')
    if not os.path.exists(co.out): 
        # make the directory thing
        os.makedirs(co.out)    

    with plt.style.context(('fivethirtyeight')):
        fig = plt.figure()  # the main animation figure
        ax = fig.add_subplot(111, aspect='equal',
                             autoscale_on=False,
                             xlim=(-20, 20),
                             ylim=(-20, 20))  # ax has all these useful attributes
        ax.grid() # Add gridlines
        no_particles = len(particles)

        pos_array = np.zeros((3, no_particles))
        vel_array = np.zeros((3, no_particles))

        if move_with == True:
            a = pos_array[0, 0]
            b = pos_array[1, 0]
            c = vel_array[0, 0]
            d = vel_array[0, 0]
        else:
            a = np.array([0,0,0])
            b = np.array([0,0,0])
            c = np.array([0,0,0])
            d = np.array([0,0,0])

        for i in range(no_particles):
            pos_array[:, i] = particles[i].x
            vel_array[:, i] = particles[i].x


        if no_particles < 10:  # Deal with multi body animate test case
            #  ax.quiver(pos_array[0, :]-a, pos_array[1, :]-b,
            #            vel_array[0, :]-c, vel_array[1, :]-d, 'g')
            ax.plot(pos_array[0, :]-a, pos_array[1, :]-b,
                    'ko', markersize=1, label='Test Masses')
            ax.plot(pos_array[0, 0]-a, pos_array[1, 0]-b, 'ko', markersize=1, label='Galaxy A')
            ax.plot(pos_array[0, -1]-a, pos_array[1, -1]-b, 'ko', markersize=1, label='Galaxy B')
        else:
            ax.plot(pos_array[0, :]-a, pos_array[1, :]-b, 'ko', markersize=0.5, label='Test Masses')
            ax.plot(pos_array[0, 0]-a, pos_array[1, 0]-b, 'g+', markersize=1, label='Galaxy A')
            #  ax.quiver(pos_array[0, -1]-a, pos_array[1, -1]-b,
            #          vel_array[0, -1]-c, vel_array[1, -1]-d)
            ax.plot(pos_array[0, -1]-a, pos_array[1, -1]-b,
                      'g+', markersize=1, label='Galaxy B')
        if no_particles > 10:
            if co.halo == True:
                if move_with == True:
                    print('I want to plot a Halo.')
        print('May have plotted Halo if needed')

        if move_with:   # In the non inertial coordinate case
            plt.xlabel(r'$X^{\prime}$')
            plt.ylabel(r'$Y^{\prime}$')
        else:
            plt.xlabel(r'$X$')
            plt.ylabel(r'$Y$')

    plt.tight_layout()
    plt.savefig(co.out + '/' + co.name + '_'+str(Time)+'t_Snap.pdf')
    plt.clf()


