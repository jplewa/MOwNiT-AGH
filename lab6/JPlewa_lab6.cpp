#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <memory>
#include <fstream>
#include <functional>
#include <cstdlib>
#include <cstdio>
#include <random>
#include <fstream>
 
namespace integration
{
	class IIntegration
	{
	    public:
        IIntegration(std::function<double(double)> function, int N) : function(function), N(N) {}
		virtual ~IIntegration() = default;

	    public:
		virtual double Integrate(double start, double end) = 0;

        double getXk(int i, double start, double end){
            return (start + (i * (end - start) / N));
        }

        protected:
        int N;
        std::function<double(double)> function;

	};

    class RectangleMethodIntegration : public IIntegration
    {
        public:
	    RectangleMethodIntegration(std::function<double(double)> function, int N) : IIntegration(function, N) {}
    
	    public:
		double Integrate(double start, double end) override
		{
            double dx = (end - start) / N;
            double rectangleSum = 0;

            for (int i = 1; i <= N; i++)
            {
                rectangleSum += function(getXk(i, start, end));
            }
            
            return (dx * rectangleSum);
        }
	};

    class TrapezoidalMethodIntegration : public IIntegration
    {
        public:
	    TrapezoidalMethodIntegration(std::function<double(double)> function, int N) : IIntegration(function, N) {}

	    public:
		double Integrate(double start, double end) override
		{
            double dx = (end - start) / N;
            double trapezoidSum = 0;
            std::vector<double> fx;

            for (int i = 0; i <= N; i++)
            {
                fx.push_back(function(getXk(i, start, end)));
            }

            for (int i = 1; i <= N; i++)
            {
                trapezoidSum += (fx[i] + fx[i - 1]);
            }

            return trapezoidSum * dx / 2;
        }
	};

    class SimpsonMethodIntegration : public IIntegration
    {
        public:
	    SimpsonMethodIntegration(std::function<double(double)> function, int N) : IIntegration(function, N) {}

	    public:
		double Integrate(double start, double end) override
		{
            double dx = (end - start) / N;
            
            double simpsonSum = 0;
            double xLeft;
            double xRight = start;

            for (int i = 1; i <= N; i++)
            {
                xLeft = xRight;
                xRight = getXk(i, start, end);
                simpsonSum += function(xLeft) + function(xRight) + 4 * function((xRight + xLeft) / 2.0);
            }

            return simpsonSum * dx / 6;
        }
	};

    class MonteCarloIntegration : public IIntegration
    {
        public:
	    MonteCarloIntegration(std::function<double(double)> function, int N) : IIntegration(function, N) {}

	    public:
		double Integrate(double start, double end) override
		{
            std::default_random_engine gen;
            std::uniform_real_distribution<> distr; 
            
            double dx = end - start;
            
            double monteCarloSum = 0;

            for (int i = 0; i < N; i++)
            {
                monteCarloSum += function(start + distr(gen) * dx);
            }

            return monteCarloSum * dx / N;
        }
	};
}

double ex1a(std::function<double(double)> fun, int N, double start, double end)
{
    std::unique_ptr<integration::RectangleMethodIntegration> rectangleMethodIntegration = 
					std::make_unique<integration::RectangleMethodIntegration>(integration::RectangleMethodIntegration(fun, N));
	return rectangleMethodIntegration -> Integrate(start, end);
}

double ex1b(std::function<double(double)> fun, int N, double start, double end)
{
    std::unique_ptr<integration::TrapezoidalMethodIntegration> trapezoidalMethodIntegration = 
					std::make_unique<integration::TrapezoidalMethodIntegration>(integration::TrapezoidalMethodIntegration(fun, N));
	return trapezoidalMethodIntegration -> Integrate(start, end);
}

double ex1c(std::function<double(double)> fun, int N, double start, double end)
{
    std::unique_ptr<integration::SimpsonMethodIntegration> simpsonMethodIntegration = 
					std::make_unique<integration::SimpsonMethodIntegration>(integration::SimpsonMethodIntegration(fun, N));
	return simpsonMethodIntegration -> Integrate(start, end);
}

double ex4(std::function<double(double)> fun, int N, double start, double end)
{
    std::unique_ptr<integration::MonteCarloIntegration> monteCarloIntegration = 
					std::make_unique<integration::MonteCarloIntegration>(integration::MonteCarloIntegration(fun, N));
	return monteCarloIntegration -> Integrate(start, end);
}

void comparison()
{
    std::ofstream comparison;
    comparison.open("csv/comparison.csv");
    std::function<double(double)> f1 = [](double x) {return x * x + 2 * x; };
    double start = -5;
    double end = 5;
    for (int i = 0; i < 7; i++)
    {
        int N = (int) pow(2, i);;
        comparison << "f1,rectangle method," << N << "," << ex1a(f1, N, start, end) << std::endl;
        comparison << "f1,trapezoidal method," << N << "," << ex1b(f1, N, start, end) << std::endl;
        comparison << "f1,Simpson method," << N << "," << ex1c(f1, N, start, end) << std::endl;
        comparison << "f1,Monte Carlo method," << N << "," << ex4(f1, N, start, end) << std::endl;
    }
    
    std::function<double(double)> f2 = [](double x) {return 7 + x * x * x + 13 * x; };
    start = -5;
    end = 5;
    for (int i = 0; i < 7; i++)
    {
        int N = (int) pow(2, i);;
        comparison << "f2,rectangle method," << N << "," << ex1a(f2, N, start, end) << std::endl;
        comparison << "f2,trapezoidal method," << N << "," << ex1b(f2, N, start, end) << std::endl;
        comparison << "f2,Simpson method," << N << "," << ex1c(f2, N, start, end) << std::endl;
        comparison << "f2,Monte Carlo method," << N << "," << ex4(f2, N, start, end) << std::endl;
    }

    std::function<double(double)> f3 = [](double x) {return 5 + x * x * x * x * x - 13 * pow(exp(1), x); };
    start = -5;
    end = 5;
    for (int i = 0; i < 10; i++)
    {
        int N = (int) pow(2, i);;
        comparison << "f3,rectangle method," << N << "," << ex1a(f3, N, start, end) << std::endl;
        comparison << "f3,trapezoidal method," << N << "," << ex1b(f3, N, start, end) << std::endl;
        comparison << "f3,Simpson method," << N << "," << ex1c(f3, N, start, end) << std::endl;
        comparison << "f3,Monte Carlo method," << N << "," << ex4(f3, N, start, end) << std::endl;
    }

    std::function<double(double)> f4 = [](double x) {return x * log(x) + pow(x, 3) - exp(x / 2) - pow(x, 2); };
    start = 0.5;
    end = 1;
    for (int i = 0; i < 10; i++)
    {
        int N = (int) pow(2, i);;
        comparison << "f4,rectangle method," << N << "," << ex1a(f4, N, start, end) << std::endl;
        comparison << "f4,trapezoidal method," << N << "," << ex1b(f4, N, start, end) << std::endl;
        comparison << "f4,Simpson method," << N << "," << ex1c(f4, N, start, end) << std::endl;
        comparison << "f4,Monte Carlo method," << N << "," << ex4(f4, N, start, end) << std::endl;
    }

    comparison.close();
}

double estimatePi(int N)
{
    std::default_random_engine gen;
    std::uniform_real_distribution<> distr; 
    int hits = 0; 
    for (int i = 0; i < N; i++) { 
        float x = distr(gen); 
        float y = distr(gen); 
        float l = sqrt(x * x + y * y); 
        if (l <= 1) 
        {
            hits++; 
        }
    } 
    return ((double) 4 * hits) / (double) N;
}

void ex3()
{
    std::ofstream ex3;
    ex3.open("csv/ex3.csv");
    for (int i = 1; i <= 1000000000; i *= 10)
    {
        ex3 << i << ',' << estimatePi(i) << std::endl;
    }
    ex3.close();
}

int main(void)
{
    ex3();
    comparison();
    return 0; 
}
