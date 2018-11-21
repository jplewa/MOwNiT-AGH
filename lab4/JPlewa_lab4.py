import numpy as np
import matplotlib.pyplot as plt
import csv
from time import time
from scipy import fftpack

def ex3():
    N = []
    dft = []
    fft = []
    _, ax = plt.subplots(figsize=(10,6))
    with open('ex3.csv', 'rb') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            N.append(int(row[0]))
            dft.append(float(row[1]))
            fft.append(float(row[2]))
    ax.set_xscale('log', basex=2)    
    ax.set_yscale('log', basey=2)
    plt.scatter(N, fft, label='FFT')
    plt.scatter(N, dft, label='DFT')
    plt.title("Speed comparison of DFT and FFT")
    ax.set_xlabel('vector size')
    ax.set_ylabel('time of execution [ms]')
    ax.legend(loc=9)
    plt.savefig("ex3.png")

def ex4():
    raw_data_y = []
    raw_data_x = []
    with open('ex4_in.csv', 'rb') as csvfile:    
        for row in csv.reader(csvfile, delimiter=','):
            if row[0] == "padding":
                break
            raw_data_x.append(float(row[0]))
            raw_data_y.append(float(row[1]))
    _, ax = plt.subplots(figsize=(10,6))
    plt.plot(raw_data_x, raw_data_y)
    plt.title("Average number of sunspots throughout the years")
    ax.set_xlabel('year')
    ax.set_ylabel('# of sunspots')
    plt.savefig("ex4_data.png")    
    re = []
    im = []
    with open('ex4_out.csv', 'rb') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            re.append(float(row[0]))
            im.append(float(row[1]))
    mag = [(r ** 2 + i ** 2) ** 0.5 for r, i in zip(re, im)]
    f_s = 1
    freqs = fftpack.fftfreq(len(mag)) * f_s
    _, ax = plt.subplots(figsize=(10,6))
    ax.stem(freqs[:len(mag)/4], mag[:len(mag)/4])
    plt.title("Signal analysis")
    ax.set_xlabel('frequency in cycles/year')
    ax.set_ylabel('frequency domain magnitude')
    ax.set_ylim(0, 7000)
    plt.savefig("ex4.png")

def ex5():
    N = []
    fft = []
    fftw = []
    _, ax = plt.subplots(figsize=(10,6))
    with open('ex5.csv', 'rb') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            N.append(int(row[0]))
            fft.append(float(row[1]))
            fftw.append(float(row[2]))
    ax.set_xscale('log', basex=2)    
    ax.set_yscale('log', basey=2)
    plt.scatter(N, fft, label='naive')
    plt.scatter(N, fftw, label='FFTW')
    plt.title("Speed comparison of naive FFT and FFTW-implementation")
    ax.set_xlabel('vector size')
    ax.set_ylabel('time of execution [ms]')
    ax.legend(loc=9)
    plt.savefig("ex5.png")

ex3()
ex4()
ex5()
