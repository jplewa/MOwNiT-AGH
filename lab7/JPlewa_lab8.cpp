#include <iostream>
#include <functional>
#include <vector>
#include <memory>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>

using namespace boost::numeric::odeint;
 
#define FIXED_POINT 100

typedef boost::array<double, 3> stateType;

namespace lorenz
{
    class LorenzSystem
    {
        public:
        LorenzSystem(double rho, double sigma, double beta) : rho(rho), sigma(sigma), beta(beta) {}

        std::vector<double> getDerivatives(double x, double y, double z)
        {
            std::vector<double> derivatives( 
                {
                    sigma * (y - x), 
                    x * (rho - z) - y, 
                    x * y - beta * z
                } 
            );
            return derivatives;
        }

        protected:
        double rho;
        double sigma;
        double beta;
    };

    class IDifferentialEquation
	{
	    public:
        IDifferentialEquation(LorenzSystem lorenzSystem, std::vector<double> initialState) : lorenzSystem(lorenzSystem), initialState(initialState) {}
		virtual ~IDifferentialEquation() = default;

	    public:
		virtual std::vector<std::vector<double>> solve(double t0, double h, int n) = 0;

        protected:
        LorenzSystem lorenzSystem;
        std::vector<double> initialState;
	};

    class Euler : IDifferentialEquation
    {
        public:
        Euler(LorenzSystem lorenzSystem, std::vector<double> initialState) : IDifferentialEquation(lorenzSystem, initialState) {}
        
        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::vector<std::vector<double>> results;
            
            double x = initialState[0];
            double y = initialState[1];
            double z = initialState[2];
            
            double t = t0;
            
            for (int j = 0; j < n; j++)
            {
                std::vector<double> derivatives = lorenzSystem.getDerivatives(x, y, z);
                x = x + h * derivatives[0];
                y = y + h * derivatives[1];
                z = z + h * derivatives[2];
                t += h;
                results.push_back({x, y, z, t});
            }
            return results;
        }
    };

    class BackwardEuler : IDifferentialEquation
    {
        public:
        BackwardEuler(LorenzSystem lorenzSystem, std::vector<double> initialState) : IDifferentialEquation(lorenzSystem, initialState) {}

        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::vector<std::vector<double>> results;

            double x = initialState[0];
            double y = initialState[1];
            double z = initialState[2];

            double t = t0;
    
            for (int j = 0; j < n; j++)
            {
                t += h;
                std::vector<double> derivatives = ApproximateDerivatives(x, y, z, h);
                x = derivatives[0];
                y = derivatives[1];
                z = derivatives[2];
                results.push_back({x, y, z, t});
            }
            return results;
        }

        private:
        std::vector<double> ApproximateDerivatives(double x, double y, double z, double h)
        {
            std::vector<double> derivatives = lorenzSystem.getDerivatives(x, y, z);
            double x0 = x + h * derivatives[0];
            double y0 = y + h * derivatives[1];
            double z0 = z + h * derivatives[2];

            derivatives = lorenzSystem.getDerivatives(x0, y0, z0);
            double x1 = x + h * derivatives[0];
            double y1 = y + h * derivatives[1];
            double z1 = z + h * derivatives[2];  

            for (int i = 0; i < FIXED_POINT - 2; i++)
            {
                derivatives = lorenzSystem.getDerivatives(x1, y1, z1);

                x0 = x1;
                y0 = y1;
                z0 = z1;

                x1 = x + h * derivatives[0];
                y1 = y + h * derivatives[1];
                z1 = z + h * derivatives[2];  
            }
            return {x1, y1, z1};
        }
    };

    class RungeKutta2 : IDifferentialEquation
    {
        public:
        RungeKutta2(LorenzSystem lorenzSystem, std::vector<double> initialState) : IDifferentialEquation(lorenzSystem, initialState) {}

        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::vector<std::vector<double>> results;

            double x = initialState[0];
            double y = initialState[1];
            double z = initialState[2];

            double t = t0;

            for (int j = 0; j < n; j++)
            {
                std::vector<double> derivatives = lorenzSystem.getDerivatives(x, y, z);
                derivatives = lorenzSystem.getDerivatives(x + h * derivatives[0] / 2, y + h * derivatives[1] / 2, z + h * derivatives[2] / 2);

                x = x + h * derivatives[0];
                y = y + h * derivatives[1];
                z = z + h * derivatives[2];
                t += h;
                results.push_back({x, y, z, t});
            }
            return results;
        }
    
    };

    class RungeKutta4 : IDifferentialEquation
    {
        public:
        RungeKutta4(LorenzSystem lorenzSystem, std::vector<double> initialState) : IDifferentialEquation(lorenzSystem, initialState) {}

        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::vector<std::vector<double>> results;

            double x = initialState[0];
            double y = initialState[1];
            double z = initialState[2];

            double t = t0;

            for (int j = 0; j < n; j++)
            {
                double k_x, k_y, k_z, sum_x, sum_y, sum_z;

                // k1
                std::vector<double> derivatives = lorenzSystem.getDerivatives(x, y, z);
                k_x = h * derivatives[0];
                k_y = h * derivatives[1];
                k_z = h * derivatives[2];
                sum_x = k_x;
                sum_y = k_y;
                sum_z = k_z;

                // k2
                derivatives = lorenzSystem.getDerivatives(x + k_x  / 2, y + k_y / 2, z + k_z / 2);
                k_x = h * derivatives[0];
                k_y = h * derivatives[1];
                k_z = h * derivatives[2];
                sum_x += 2 * k_x;
                sum_y += 2 * k_y;
                sum_z += 2 * k_z;

                // k3
                derivatives = lorenzSystem.getDerivatives(x + k_x  / 2, y + k_y / 2, z + k_z / 2);
                k_x = h * derivatives[0];
                k_y = h * derivatives[1];
                k_z = h * derivatives[2];
                sum_x += 2 * k_x;
                sum_y += 2 * k_y;
                sum_z += 2 * k_z;

                // k4
                derivatives = lorenzSystem.getDerivatives(x + k_x, y + k_y, z + k_z);
                k_x = h * derivatives[0];
                k_z = h * derivatives[2];
                k_y = h * derivatives[1];
                sum_x += k_x;
                sum_y += k_y;
                sum_z += k_z;

                x = x + sum_x / 6;
                y = y + sum_y / 6;
                z = z + sum_z / 6;
                
                t += h;
                results.push_back({x, y, z, t});
            }
            return results;
        }
    };

    class BoostSolver : IDifferentialEquation
    {
        public:
        BoostSolver(LorenzSystem lorenzSystem, std::vector<double> initialState) : IDifferentialEquation(lorenzSystem, initialState) {}

        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::function<void(stateType, stateType&, double)> der = [&](stateType x, stateType &dxdt, double t)
            {
                std::vector<double> derivatives = this -> lorenzSystem.getDerivatives(x[0], x[1], x[2]);
                dxdt[0] = derivatives[0];
                dxdt[1] = derivatives[1];
                dxdt[2] = derivatives[2];
            };

            std::function<void(stateType, double)> save = [&](stateType x, double t)
            {
                appendSolution(x, t);
            };

            stateType state = {initialState[0], initialState[1], initialState[2]};

            integrate(der, state, t0, (double) n * h, h, save);

            return this -> solution;
        }

        void appendSolution(stateType x, double t)
        {
            solution.push_back({x[0], x[1], x[2], t}); 
        }

        private:
        std::vector<std::vector<double>> solution;
    };
}

void save_solution(std::vector<std::vector<double>> result, std::ofstream &file)
{
    for (int i = 0; i < result.size(); i++)
    {
        for (int j = 0; j < result[i].size(); j++)
        {
            file << result[i][j] << ",";
        }
        file << std::endl;
    }
}

void ex2a()
{
    std::ofstream file;
    file.open("csv/euler.csv");
    lorenz::LorenzSystem lorenz(28, 10, 8.0 / 3.0);
    std::vector<double> state({1, 1, 1});
    std::unique_ptr<lorenz::Euler> solver = std::make_unique<lorenz::Euler>(lorenz::Euler(lorenz, state));    
    std::vector<std::vector<double>> result = solver -> solve(0, 0.01, 4000);
    save_solution(result, file);
    file.close();
}

void ex2b()
{
    std::ofstream file;
    file.open("csv/backward_euler.csv");
    lorenz::LorenzSystem lorenz(28, 10, 8.0 / 3.0);
    std::vector<double> state({1, 1, 1});
    std::unique_ptr<lorenz::BackwardEuler> solver = std::make_unique<lorenz::BackwardEuler>(lorenz::BackwardEuler(lorenz, state));
    std::vector<std::vector<double>> result = solver -> solve(0, 0.01, 4000);
    save_solution(result, file);
    file.close();
}

void ex2d()
{
    std::ofstream file;
    file.open("csv/rk2.csv");
    lorenz::LorenzSystem lorenz(28, 10, 8.0 / 3.0);
    std::vector<double> state({1, 1, 1});
    std::unique_ptr<lorenz::RungeKutta2> solver = std::make_unique<lorenz::RungeKutta2>(lorenz::RungeKutta2(lorenz, state));
    std::vector<std::vector<double>> result = solver -> solve(0, 0.01, 4000);
    save_solution(result, file);
    file.close();
}

void ex2e()
{
    std::ofstream file;
    file.open("csv/rk4.csv");
    lorenz::LorenzSystem lorenz(28, 10, 8.0 / 3.0);
    std::vector<double> state({1, 1, 1});
    std::unique_ptr<lorenz::RungeKutta4> solver = std::make_unique<lorenz::RungeKutta4>(lorenz::RungeKutta4(lorenz, state));
    std::vector<std::vector<double>> result = solver -> solve(0, 0.01, 4000);
    save_solution(result, file);
    file.close();
}

void ex4()
{
    std::ofstream file;
    file.open("csv/boost.csv");
    lorenz::LorenzSystem lorenz(28, 10, 8.0 / 3.0);
    std::vector<double> state({1, 1, 1});
    std::unique_ptr<lorenz::BoostSolver> solver = std::make_unique<lorenz::BoostSolver>(lorenz::BoostSolver(lorenz, state));
    std::vector<std::vector<double>> result = solver -> solve(0, 0.01, 4000);
    save_solution(result, file);
    file.close();
}

int main(void)
{
    ex2a();
    ex2b();
    ex2d();
    ex2e();
    ex4();
    return 0;
}