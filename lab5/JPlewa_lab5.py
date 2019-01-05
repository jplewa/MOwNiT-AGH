import csv
import collections
import matplotlib.pyplot as plt


def ex5():
    jacobi_iter = collections.defaultdict(list)
    jacobi_conv = collections.defaultdict(list)

    gauss_seidel_iter = collections.defaultdict(list)
    gauss_seidel_conv = collections.defaultdict(list)

    sor_iter = collections.defaultdict(list)
    sor_conv = collections.defaultdict(list)

    with open('csv/ex5_jacobi.csv', 'rt') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            jacobi_iter[int(row[0])].append(int(row[1]))
            jacobi_conv[int(row[0])].append(float(row[2]))

    with open('csv/ex5_gauss_seidel.csv', 'rt') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            gauss_seidel_iter[int(row[0])].append(int(row[1]))
            gauss_seidel_conv[int(row[0])].append(float(row[2]))

    with open('csv/ex5_sor.csv', 'rt') as csvfile:
        for row in csv.reader(csvfile, delimiter=','):
            sor_iter[int(row[0])].append(int(row[1]))
            sor_conv[int(row[0])].append(float(row[2]))

    for i in range(len(jacobi_iter.keys())):
        _, ax = plt.subplots(figsize=(10, 6))
        plt.plot(jacobi_iter[i], jacobi_conv[i], label='jacobi')
        plt.plot(sor_iter[i], sor_conv[i], label='sor')
        plt.plot(gauss_seidel_iter[i],
                 gauss_seidel_conv[i], label='gauss_seidel')
        ax.legend(loc=9)
        plt.title('Iterative method convergence comparison')
        ax.set_xlabel('iteration')
        ax.set_ylabel(r'$\max_{i}(|x^{(i)} - x^{(i-1)}|)$')
        ax.set_yscale('log', basey=10)
        plt.savefig(f'plot/ex5_{i}.png')


ex5()
