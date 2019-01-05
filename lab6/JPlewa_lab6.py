import matplotlib.pyplot as plt
from collections import defaultdict
import csv
import math


def comparison():
    expected = {}
    with open('JPlewa_comparison_expected.csv', 'rt') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            expected[row[0]] = float(row[1])
    #x = [2 ** i for i in range(0, 25)]
    xs = defaultdict(lambda: defaultdict(list))
    ys = defaultdict(lambda: defaultdict(list))
    with open('comparison.csv', 'rt') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            xs[row[0]][row[1]].append(float(row[2]))
            ys[row[0]][row[1]].append(float(row[3]))
    for fun, methods in ys.items():
        print(fun)
        _, ax = plt.subplots(figsize=(10,6))
        for method, values in methods.items():
            plt.scatter(xs[fun][method], values, label=method)
        plt.axhline(expected[fun], label='exact', color='gray')
        ax.set_xscale('log', basex=2)    
        plt.title('Convergence of integration methods')
        ax.set_xlabel('N')
        ax.set_ylabel('result')
        ax.legend()
        plt.savefig(f'{fun}.png')


def ex3():
    xs = []
    ys = []
    with open('ex3.csv', 'rt') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            xs.append(float(row[0]))
            ys.append(float(row[1]))
    _, ax = plt.subplots(figsize=(10,6))
    plt.scatter(xs, ys, label='estimation')
    plt.axhline(math.pi, label='exact', color='gray')
    ax.set_xscale('log', basex=10)    
    plt.title('Convergence of Monte Carlo PI estimation')
    ax.set_xlabel('N')
    ax.set_ylabel('pi estimation')
    ax.legend()
    plt.savefig(f'ex3.png')

comparison()
ex3()