import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv
import math
import os
import sys


LORENZ_FILES = ['euler.csv', 'backward_euler.csv',
                'rk2.csv', 'rk4.csv', 'boost.csv']


def plot_lorenz(csv_name):
    x = []
    y = []
    z = []
    with open(os.path.join('csv', csv_name), 'rt') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            x.append(float(row[0]))
            y.append(float(row[1]))
            z.append(float(row[2]))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(x, y, z, label=csv_name.split('.')[0])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.title('Lorenz system solution', fontsize=15, y=1.08)
    mpl.rcParams['legend.fontsize'] = 10
    ax.legend()
    plt.savefig(os.path.join('plot', csv_name.replace('.csv', '.png')), quality=100, dpi=500)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        for file in LORENZ_FILES:
            plot_lorenz(file)
    elif sys.argv[1].isdigit():
        print(int(sys.argv[1]))
