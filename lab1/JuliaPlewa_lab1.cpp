#include <stdio.h>
#include <math.h>
#include <gsl/gsl_ieee_utils.h>
#include <iostream>
#include <fstream>

#define EX2_LEN 9999    // exercise 2 sequence length
#define EX3_LEN 2000    // exercise 3 sequence length

void print_to_files(float a, double b, std::ofstream *file_a, std::ofstream *file_b) {
    *file_a << a << "\n";
    *file_b << b << "\n";
}

void ex2() {
    std::ofstream file_a;   // open result files and set precision
    std::ofstream file_b;
    file_a.precision(6);
    file_b.precision(15);
    file_a.setf(std::ios::scientific);
    file_b.setf(std::ios::scientific);
    file_a.open("ex2_a.txt");
    file_b.open("ex2_b.txt");

    float a;                // initial sequence values
    float prev_a = 30;
    float prev_prev_a = 15;
    double b;
    double prev_b = 30;
    double prev_prev_b = 15;

    // save first elements to their respective files
    print_to_files(prev_prev_a, prev_prev_b, &file_a, &file_b);
    print_to_files(prev_a, prev_b, &file_a, &file_b);

    // loop through sequence elements and save values into files
    for (int i = 0; i < EX2_LEN; i++) {
        a = prev_a * (float) 6 / (float) 7 - prev_prev_a * (float) 1 / (float) 7;
        b = prev_b * (double) 6 / (double) 7 - prev_prev_b * (double) 1 / (double) 7;
        prev_prev_a = prev_a;
        prev_prev_b = prev_b;
        prev_a = a;
        prev_b = b;
        print_to_files(a, b, &file_a, &file_b);
    }
    file_a.close();     // close result files
    file_b.close();
}

void ex3() {
    // open result files
    FILE *file_a = fopen("ex3_a.txt", "w+");
    FILE *file_b = fopen("ex3_b.txt", "w+");

    float a = 1.1;      // set initial sequence values
    double b = 1.1;

    // loop through the sequence and record the elemts
    for (int i = 0; i < EX3_LEN; i++) {
        gsl_ieee_fprintf_float(file_a, &a);
        gsl_ieee_fprintf_double(file_b, &b);
        fprintf(file_a, "\n");
        fprintf(file_b, "\n");
        a /= ((float) 2);
        b /= ((double) 2);
    }
    fclose(file_a);     // close result files
    fclose(file_b);
}

void ex4() {
    std::cout.precision(15);    // set output precision
    double values[] = {0.00000000000001,  0.00000001826, 1.0, 10.0};
    // loop through the sample x values and calculate f(x) and g(x) separately
    // printing out the results
    for (int i = 0; i < 4; i++) {
        double x = values[i];
        double fx = sqrt(x * x + 1) - 1;
        double gx = (x * x) / (sqrt(x * x + 1) + 1);

        std::cout << "x = " << x << std::endl << "\tf(x) = " << fx << "\t";
        gsl_ieee_fprintf_double(stdout, &fx);
        std::cout << std::endl << "\tg(x) = " << gx << "\t";
        gsl_ieee_fprintf_double(stdout, &gx);
        std::cout << std::endl << std::endl;
    }
}

int main() {
    std::cout.setf(std::ios::scientific);
    ex2();
    ex3();
    ex4();
    return 0;
}
