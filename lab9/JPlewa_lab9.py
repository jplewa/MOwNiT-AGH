import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import math
import os
import sys

ODE_FILES = ['euler.csv', 'backward_euler.csv',
             'rk2.csv', 'rk4.csv', 'boost.csv']

PREDATOR_PREY = 'predator_prey'
MATHEMATICAL_PENDULUM = 'mathematical_pendulum'
SPRING = 'spring'
DECAY = 'decay'

PLOT_DIR = 'plot'
CSV_DIR = 'csv'


def predator_prey():
    for ode_file in ODE_FILES:
        t = []
        prey = []
        predators = []
        with open(os.path.join(CSV_DIR, os.path.join(PREDATOR_PREY, ode_file)), 'rt') as csvfile:
            for row in csv.reader(csvfile, delimiter=','):
                predators.append(float(row[0]))
                prey.append(float(row[1]))
                t.append(float(row[2]))
        _, ax = plt.subplots(figsize=(10, 6))
        plt.plot(t, predators, label='predators')
        plt.plot(t, prey, label='prey')
        ax.legend(loc=9)
        plt.title(
            ('Prey-predators system solution: ' + ode_file.split('.')[0]))
        ax.set_xlabel('time')
        plt.savefig(os.path.join(PLOT_DIR, os.path.join(
            PREDATOR_PREY, ode_file.replace('.csv', '.png'))), quality=100, dpi=500)

        _, ax = plt.subplots(figsize=(10, 6))
        plt.plot(prey, predators)
        plt.title(
            'Phase space plot of prey-predators system solution: ' + ode_file.split('.')[0])
        ax.set_xlabel('prey')
        ax.set_ylabel('predators')
        plt.savefig(os.path.join(PLOT_DIR, os.path.join(
            PREDATOR_PREY, 'phase_' + ode_file.replace('.csv', '.png'))), quality=100, dpi=400)


def mathematical_pendulum():
    for ode_file in ODE_FILES:
        t = []
        omega = []
        theta = []
        with open(os.path.join(CSV_DIR, os.path.join(MATHEMATICAL_PENDULUM, ode_file)), 'rt') as csvfile:
            for row in csv.reader(csvfile, delimiter=','):
                theta.append(math.degrees(float(row[0])))
                omega.append(math.degrees(float(row[1])))
                t.append(float(row[2]))
        _, ax = plt.subplots(figsize=(10, 6))
        plt.plot(t, theta, label='theta')
        plt.plot(t, omega, label='omega')
        ax.legend(loc=9)
        plt.title(f'Mathematical pendulum solution: {ode_file.split(".")[0]}')
        ax.set_xlabel('time')
        ax.set_ylabel('degrees')
        plt.savefig(os.path.join(PLOT_DIR, os.path.join(
            MATHEMATICAL_PENDULUM, ode_file.replace('.csv', '.png'))), quality=100, dpi=400)


def spring():
    for ode_file in ODE_FILES:
        y1 = []
        y2 = []
        t = []
        with open(os.path.join(CSV_DIR, os.path.join(SPRING, ode_file)), 'rt') as csvfile:
            for row in csv.reader(csvfile, delimiter=','):
                y1.append(float(row[0]))
                y2.append(float(row[1]))
                t.append(float(row[2]))
        _, ax = plt.subplots(figsize=(10, 6))
        plt.plot(t, y2, label='x')
        ax.legend(loc=9)
        plt.title('Spring: ' + ode_file.split('.')[0])
        ax.set_xlabel('time')
        plt.savefig(os.path.join(PLOT_DIR, os.path.join(
            SPRING, ode_file.replace('.csv', '.png'))), quality=100, dpi=400)


def decay():
    for ode_file in ODE_FILES:
        u = []
        t = []
        with open(os.path.join(CSV_DIR, os.path.join(DECAY, ode_file)), 'rt') as csvfile:
            for row in csv.reader(csvfile, delimiter=','):
                u.append(float(row[0]))
                t.append(float(row[1]))
        _, ax = plt.subplots(figsize=(10, 6))
        plt.plot(t, u, label='u')
        plt.plot(t, [math.e ** (-25 * x) for x in t], label='u_analytical')
        ax.legend(loc=9)
        plt.title('Radioactive decay: ' + ode_file.split('.')[0])
        ax.set_xlabel('time')
        plt.savefig(os.path.join(PLOT_DIR, os.path.join(
            DECAY, ode_file.replace('.csv', '.png'))), quality=100, dpi=400)


if __name__ == "__main__":
    # predator_prey()
    # mathematical_pendulum()
    spring()
    # decay()
