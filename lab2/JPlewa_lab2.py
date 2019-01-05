import matplotlib.pyplot as plt


def ex2():
    with open("csv/cubic.csv", "r") as cubic_file:
        xs = []
        cubic = []
        for line in cubic_file.readlines():
            xs.append(float(line.split(',')[0]))
            cubic.append(float(line.split(',')[1]))

    with open("csv/lagrange.csv", "r") as lagrange_file:
        lagrange = []
        for line in lagrange_file.readlines():
            lagrange.append(float(line.split(',')[1]))

    _, ax = plt.subplots(figsize=(12, 5))
    ax.plot(xs, lagrange, label='Lagrange')
    ax.plot(xs, cubic, label='cubic Hermite spline')
    ax.legend(loc='upper center', shadow=True, fontsize='x-large')
    plt.ylabel('f(x)')
    plt.xlabel('x')
    plt.savefig('plot/ex2.png')


xs = [108, 19, 13, 124, 40, 57, 23, 14, 45, 10,
      5, 48, 11, 23, 7, 2, 24, 6, 3, 23, 6, 9,
      9, 3, 29, 7, 4, 20, 7, 4, 0, 25, 6, 5, 22,
      11, 61, 12, 4, 16, 13, 60, 41, 37, 55, 41,
      11, 27, 8, 3, 17, 13, 13, 15, 8, 29, 30,
      24, 9, 31, 14, 53, 26]

ys = [392.5, 46.2, 15.7, 422.2, 119.4, 170.9, 56.9,
      77.5, 214, 65.3, 20.9, 248.1, 23.5, 39.6, 48.8,
      6.6, 134.9, 50.9, 4.4, 113, 14.8, 48.7, 52.1,
      13.2, 103.9, 77.5, 11.8, 98.1, 27.9, 38.1,
      0, 69.2, 14.6, 40.3, 161.5, 57.2, 217.6,
      58.1, 12.6, 59.6, 89.9, 202.4, 181.3,
      152.8, 162.8, 73.4, 21.3, 92.6, 76.1,
      39.9, 142.1, 93, 31.9, 32.1, 55.6, 133.3,
      194.5, 137.9, 87.4, 209.8, 95.5, 244.6, 187.5]


def ex5():
    _, _ = plt.subplots(figsize=(12, 5))
    xs_ys_sorted = list(zip(*sorted((zip(xs, ys)))))
    xs_sorted = xs_ys_sorted[0]
    ys_sorted = xs_ys_sorted[1]
    ys_predicted = [3.41382 * x + 19.9945 for x in xs_ys_sorted[0]]
    plt.scatter(xs_sorted, ys_sorted, marker='.')
    plt.plot(xs_sorted, ys_predicted, color='orange')
    plt.ylabel('f(x)')
    plt.xlabel('x')
    plt.savefig('plot/ex5.png')


def ex6():
    predicted_ys = [388.687, 84.8571, 64.3742, 443.309, 156.547,
                    214.582, 98.5124, 67.788, 173.617, 54.1327,
                    37.0636, 183.858, 57.5465, 98.5124, 43.8912,
                    26.8221, 101.926, 40.4774, 30.2359, 98.5124,
                    40.4774, 50.7189, 50.7189, 30.2359, 118.995,
                    43.8912, 33.6498, 88.271, 43.8912, 33.6498,
                    19.9945, 105.34, 40.4774, 37.0636, 95.0986,
                    57.5465, 228.238, 60.9604, 33.6498, 74.6157,
                    64.3742, 224.824, 159.961, 146.306, 207.755,
                    159.961, 57.5465, 112.168, 47.3051, 30.2359,
                    78.0295, 64.3742, 64.3742, 71.2018, 47.3051,
                    118.995, 122.409, 101.926, 50.7189, 125.823,
                    67.788, 200.927, 108.754]
    ax, _ = plt.subplots(figsize=(12, 5))
    xs_ys_sorted = list(zip(*sorted((zip(xs, ys)))))
    xs_sorted = xs_ys_sorted[0]
    ys_sorted = xs_ys_sorted[1]
    predicted_ys_sorted = list(zip(*sorted((zip(xs, predicted_ys)))))[1]
    plt.scatter(xs_sorted, ys_sorted, marker="^", label='input y')
    plt.scatter(xs_sorted, predicted_ys_sorted, marker='*',
                color='orange', label='predicted y')
    ax.legend(loc='center right', shadow=True, fontsize='x-large')
    plt.ylabel('f(x)')
    plt.xlabel('x')
    plt.savefig('plot/ex6.png')


ex2()
ex5()
ex6()
