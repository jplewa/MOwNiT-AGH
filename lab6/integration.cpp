#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include <algorithm>
#include <memory>
#include <fstream>
#include <functional>

#define N 10

namespace integration
{
	class IIntegration
	{
	    public:
		virtual ~IIntegration() = default;

	    public:
		virtual double Integrate(double start, double end) = 0;

        double getXk(int i, double start, double end){
            return (start + (i * (end - start) / N));
        }

	};

    class RectangleMethodIntegration : public IIntegration
    {
        public:
	    RectangleMethodIntegration(std::function<double(double)> function) : function(function) {}

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
	    
        private:
        std::function<double(double)> function;
	};

    class TrapezoidalMethodIntegration : public IIntegration
    {
        public:
	    TrapezoidalMethodIntegration(std::function<double(double)> function) : function(function) {}

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
	    
        private:
        std::function<double(double)> function;
	};

    class SimpsonMethodIntegration : public IIntegration
    {
        public:
	    SimpsonMethodIntegration(std::function<double(double)> function) : function(function) {}

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
	    
        private:
        std::function<double(double)> function;
	};
}

void ex1a()
{
    std::function<double(double)> test = [](double x) { return x; };

    std::unique_ptr<integration::RectangleMethodIntegration> rectangleMethodIntegration = 
					std::make_unique<integration::RectangleMethodIntegration>(integration::RectangleMethodIntegration(test));
	std::cout << rectangleMethodIntegration -> Integrate(0, 2);
	std::cout << std::endl;
}

void ex1b()
{
    std::function<double(double)> test = [](double x) { return x; };

    std::unique_ptr<integration::TrapezoidalMethodIntegration> trapezoidalMethodIntegration = 
					std::make_unique<integration::TrapezoidalMethodIntegration>(integration::TrapezoidalMethodIntegration(test));
	std::cout << trapezoidalMethodIntegration -> Integrate(0, 2);
	std::cout << std::endl;
}

void ex1c()
{
    std::function<double(double)> test = [](double x) { return x * x + 2 * x; };

    std::unique_ptr<integration::SimpsonMethodIntegration> simpsonMethodIntegration = 
					std::make_unique<integration::SimpsonMethodIntegration>(integration::SimpsonMethodIntegration(test));
	std::cout << simpsonMethodIntegration -> Integrate(0, 1);
	std::cout << std::endl;
}

int main(void)
{
    ex1c();
    return 0;
}