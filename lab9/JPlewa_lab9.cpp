#include <math.h>
#include <stdio.h>
#include <iostream>
#include <functional>
#include <vector>
#include <memory>
#include <fstream>
#include <boost/numeric/odeint.hpp>
#include <boost/array.hpp>

using namespace boost::numeric::odeint;
 
#define FIXED_POINT 100
#define MAX_DIMENSIONS 4
#define PI 3.14159265

typedef boost::array<double, MAX_DIMENSIONS> stateType;

std::vector<std::vector<double>> leapFrogAlgorithm(std::vector<double> a, double x0, double v0, double start, double dt, int N)
{
    std::vector<std::vector<double>> results;

    std::vector<double> t(N);
    std::vector<double> x(N);
    std::vector<double> v(N);

    t[0] = start;
    x[0] = x0;
    v[0] = v0;

    results.push_back({t[0], x[0], v[0], a[0]});


    for (int i = 0; i < N - 1; i++)
    {
        t[i + 1] = t[i] + dt;
        x[i + 1] = x[i] + v[i] * dt + dt * dt * a[i] / 2;
        v[i + 1] = v[i] + dt * (a[i] + a[i + 1]) / 2;
        results.push_back({t[i + 1], x[i + 1], v[i + 1], a[i + 1]});
    }
    return results;
}

namespace ode
{
    class IDifferentialEquation
    {
        public:
        IDifferentialEquation(std::vector<double> constants) : constants(constants) {}
		virtual ~IDifferentialEquation() = default;

        public:
        virtual std::vector<double> getDerivatives(std::vector<double> variables) = 0;
        virtual int getDimensions() = 0;

        protected:
        std::vector<double> constants;
    };

    class PredatorPreyODE : public IDifferentialEquation
    {
        public:
        PredatorPreyODE(std::vector<double> constants) : IDifferentialEquation(constants) {}

        std::vector<double> getDerivatives(std::vector<double> variables)
        {
            double x = variables[0];
            double y = variables[1];

            std::vector<double> derivatives({
                getC() * getD() * x * y - getA() * x,
                getB() * y - getD() * x * y
            });

            return derivatives;
        }
        int getDimensions()
        {
            return 2;
        }

        private:
        double getA()
        {
            return constants[0];
        }
        double getB()
        {
            return constants[1];
        }
        double getC()
        {
            return constants[2];
        }
        double getD()
        {
            return constants[3];
        }
    };

    class MathematicalPendulumODE : public IDifferentialEquation
    {
        public:
        MathematicalPendulumODE(std::vector<double> constants) : IDifferentialEquation(constants) {}

        std::vector<double> getDerivatives(std::vector<double> variables)
        {
            double omega = variables[0];
            double theta = variables[1];

            std::vector<double> derivatives({
                -9.80665 * sin(theta) / getL(),
                omega
            });

            return derivatives;
        }
        int getDimensions()
        {
            return 2;
        }

        private:
        double getL()
        {
            return constants[0];
        }
    };

    class SpringODE : public IDifferentialEquation
    {
        public:
        SpringODE(std::vector<double> constants) : IDifferentialEquation(constants) {}

        std::vector<double> getDerivatives(std::vector<double> variables)
        {
            double y1 = variables[0];
            double y2 = variables[1];

            std::vector<double> derivatives({
                y2,
                (-1) * y1 * getK() / getM()
            });

            return derivatives;
        }
        int getDimensions()
        {
            return 2;
        }

        private:
        double getK()
        {
            return constants[0];
        }
        double getM()
        {
            return constants[1];
        }
    };

    class DecayODE : public IDifferentialEquation
    {
        public:
        DecayODE(std::vector<double> constants) : IDifferentialEquation(constants) {}

        std::vector<double> getDerivatives(std::vector<double> variables)
        {
            double u = variables[0];

            std::vector<double> derivatives({
                (-1) * u / getTau()
            });

            return derivatives;
        }
        int getDimensions()
        {
            return 1;
        }

        private:
        double getTau()
        {
            return constants[0];
        }
    };

    class IDifferentialEquationSolver
	{
	    public:
        IDifferentialEquationSolver(IDifferentialEquation* differentialEquation, std::vector<double> initialState) : differentialEquation(differentialEquation), initialState(initialState) {}
		virtual ~IDifferentialEquationSolver() = default;

	    public:
		virtual std::vector<std::vector<double>> solve(double t0, double h, int n) = 0;

        protected:
        IDifferentialEquation* differentialEquation;
        std::vector<double> initialState;
	};

    class Euler : IDifferentialEquationSolver
    {
        public:
        Euler(IDifferentialEquation* differentialEquation, std::vector<double> initialState) : IDifferentialEquationSolver(differentialEquation, initialState) {}
        
        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::vector<std::vector<double>> results;
            
            double t = t0;

            std::vector<double> state = initialState;
            state.push_back(t);
            
            for (int j = 0; j < n; j++)
            {
                std::vector<double> derivatives = differentialEquation -> getDerivatives(state);
                for (int i = 0; i < derivatives.size(); i++)
                {
                    state[i] = state[i] + h * derivatives[i];
                }
                t += h;
                state[state.size() - 1] = t;
                results.push_back(state);
            }
            return results;
        }
    };

    class BackwardEuler : IDifferentialEquationSolver
    {
        public:
        BackwardEuler(IDifferentialEquation* differentialEquation, std::vector<double> initialState) : IDifferentialEquationSolver(differentialEquation, initialState) {}

        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::vector<std::vector<double>> results;

            double t = t0;
            std::vector<double> state = initialState;
            state.push_back(t);
    
            for (int j = 0; j < n; j++)
            {
                t += h;
                std::vector<double> derivatives = ApproximateDerivatives(state, h);
                for (int i = 0; i < derivatives.size(); i++)
                {
                    state[i] = derivatives[i];
                }
                state[state.size() - 1] = t;
                results.push_back(state);
            }
            return results;
        }

        private:
        std::vector<double> ApproximateDerivatives(std::vector<double> state, double h)
        {
            std::vector<double> state0(state.size() - 1);
            std::vector<double> state1(state.size() - 1);

            std::vector<double> derivatives = differentialEquation -> getDerivatives(state);
            for (int i = 0; i < state0.size(); i++)
            {
                state0[i] = state[i] + h * derivatives[i];
            }

            derivatives = differentialEquation -> getDerivatives(state0);
            for (int i = 0; i < state1.size(); i++)
            {
                state1[i] = state[i] + h * derivatives[i];
            }


            for (int j = 0; j < FIXED_POINT - 2; j++)
            {
                derivatives = differentialEquation -> getDerivatives(state1);

                for (int i = 0; i < state0.size(); i++)
                {
                    state0[i] = state1[i];
                }
                for (int i = 0; i < state1.size(); i++)
                {
                    state1[i] = state[i] + h * derivatives[i];
                }
            }
            return state1;
        }
    };  

    class RungeKutta2 : IDifferentialEquationSolver
    {
        public:
        RungeKutta2(IDifferentialEquation* differentialEquation, std::vector<double> initialState) : IDifferentialEquationSolver(differentialEquation, initialState) {}

        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::vector<std::vector<double>> results;

            double t = t0;

            std::vector<double> state = initialState;
            state.push_back(t);

            for (int j = 0; j < n; j++)
            {
                std::vector<double> derivatives = differentialEquation -> getDerivatives(state);
                
                std::vector<double> state1(derivatives.size());
                for (int i = 0; i < state1.size(); i++)
                {
                    state1[i] = state[i] + h * derivatives[i] / 2;
                }
                
                derivatives = differentialEquation -> getDerivatives(state1);
                for (int i = 0; i < state1.size(); i++)
                {
                    state[i] += h * derivatives[i];
                }
                t += h;
                state[state.size() - 1] = t;
                results.push_back(state);
            }
            return results;
        }
    
    };

    class RungeKutta4 : IDifferentialEquationSolver
    {
        public:
        RungeKutta4(IDifferentialEquation* differentialEquation, std::vector<double> initialState) : IDifferentialEquationSolver(differentialEquation, initialState) {}

        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::vector<std::vector<double>> results;

            double t = t0;

            std::vector<double> state = initialState;
            state.push_back(t);

            for (int j = 0; j < n; j++)
            {
                std::vector<double> k(initialState.size());
                std::vector<double> sum(initialState.size());
                std::vector<double> variables(initialState.size());

                // k1
                std::vector<double> derivatives = differentialEquation -> getDerivatives(state);
                for (int i = 0; i < derivatives.size(); i++)
                {
                    k[i] = h * derivatives[i];
                    sum[i] = k[i];
                }

                // k2, k3
                for (int l = 0; l < 2; l++)
                {
                    for (int i = 0; i < derivatives.size(); i++)
                    {
                        variables[i] = state[i]  + k[i] / 2;
                    }
                    derivatives = differentialEquation -> getDerivatives(variables);
                    for (int i = 0; i < derivatives.size(); i++)
                    {
                        k[i] = h * derivatives[i];
                        sum[i] += 2 * k[i];
                    }
                }

                // k4
                for (int i = 0; i < derivatives.size(); i++)
                {
                    variables[i] = state[i]  + k[i];
                }
                derivatives = differentialEquation -> getDerivatives(variables);
                for (int i = 0; i < derivatives.size(); i++)
                {
                    k[i] = h * derivatives[i];
                    sum[i] += k[i];
                }

                for (int i = 0; i < derivatives.size(); i++)
                {
                    state[i] += sum[i] / 6;
                }

                t += h;
                state[state.size() - 1] = t;                
                results.push_back(state);
            }
            return results;
        }
    };

    class BoostSolver : IDifferentialEquationSolver
    {
        public:
        BoostSolver(IDifferentialEquation* differentialEquation, std::vector<double> initialState) : IDifferentialEquationSolver(differentialEquation, initialState) {}

        std::vector<std::vector<double>> solve(double t0, double h, int n)
        {
            std::function<void(stateType, stateType&, double)> der = [&](stateType x, stateType &dxdt, double t)
            {
                std::vector<double> state(differentialEquation -> getDimensions());
                for (int i = 0; i < state.size(); i++)
                {
                    state[i] = x[i];
                }
                std::vector<double> derivatives = this -> differentialEquation -> getDerivatives(state);
                for (int i = 0; i < derivatives.size(); i++)
                {
                    dxdt[i] = derivatives[i];
                }
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
            std::vector<double> result(differentialEquation -> getDimensions());
            for (int i = 0; i < result.size(); i++)
            {
                result[i] = x[i];
            }
            result.push_back(t);
            solution.push_back(result); 
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

void ex2()
{
    std::ofstream file;
    file.open("csv/predator_prey/euler.csv");
    ode::PredatorPreyODE ode({0.1, 0.02, 0.4, 0.02});
    std::vector<double> state({10, 10});
    std::unique_ptr<ode::Euler> solver = std::make_unique<ode::Euler>(ode::Euler(&ode, state));    
    std::vector<std::vector<double>> result = solver -> solve(0, 0.5, 100000);
    save_solution(result, file);
    file.close();

    file.open("csv/predator_prey/backward_euler.csv");
    std::unique_ptr<ode::BackwardEuler> solver2 = std::make_unique<ode::BackwardEuler>(ode::BackwardEuler(&ode, state));    
    std::vector<std::vector<double>> result2 = solver2 -> solve(0, 0.5, 100000);
    save_solution(result2, file);
    file.close();

    file.open("csv/predator_prey/rk2.csv");
    std::unique_ptr<ode::RungeKutta2> solver3 = std::make_unique<ode::RungeKutta2>(ode::RungeKutta2(&ode, state));    
    std::vector<std::vector<double>> result3 = solver3 -> solve(0, 0.5, 100000);
    save_solution(result3, file);
    file.close();

    file.open("csv/predator_prey/rk4.csv");
    std::unique_ptr<ode::RungeKutta4> solver4 = std::make_unique<ode::RungeKutta4>(ode::RungeKutta4(&ode, state));    
    std::vector<std::vector<double>> result4 = solver4 -> solve(0, 0.5, 100000);
    save_solution(result3, file);
    file.close();

    file.open("csv/predator_prey/boost.csv");
    std::unique_ptr<ode::BoostSolver> solver5 = std::make_unique<ode::BoostSolver>(ode::BoostSolver(&ode, state));    
    std::vector<std::vector<double>> result5 = solver5 -> solve(0, 0.5, 100000);
    save_solution(result4, file);
    file.close();
}

void ex3()
{
    std::ofstream file;
    file.open("csv/mathematical_pendulum/euler.csv");
    ode::MathematicalPendulumODE ode({1});
    std::vector<double> state({0, 3.141592653e-06});
    std::unique_ptr<ode::Euler> solver = std::make_unique<ode::Euler>(ode::Euler(&ode, state));    
    std::vector<std::vector<double>> result = solver -> solve(0, 0.1, 100);
    save_solution(result, file);
    file.close();

    file.open("csv/mathematical_pendulum/backward_euler.csv");
    std::unique_ptr<ode::BackwardEuler> solver2 = std::make_unique<ode::BackwardEuler>(ode::BackwardEuler(&ode, state));    
    std::vector<std::vector<double>> result2 = solver2 -> solve(0, 0.1, 100);
    save_solution(result2, file);
    file.close();

    file.open("csv/mathematical_pendulum/rk2.csv");
    std::unique_ptr<ode::RungeKutta2> solver3 = std::make_unique<ode::RungeKutta2>(ode::RungeKutta2(&ode, state));    
    std::vector<std::vector<double>> result3 = solver3 -> solve(0, 0.1, 100);
    save_solution(result3, file);
    file.close();

    file.open("csv/mathematical_pendulum/rk4.csv");
    std::unique_ptr<ode::RungeKutta4> solver4 = std::make_unique<ode::RungeKutta4>(ode::RungeKutta4(&ode, state));    
    std::vector<std::vector<double>> result4 = solver4 -> solve(0, 0.1, 100);
    save_solution(result3, file);
    file.close();

    file.open("csv/mathematical_pendulum/boost.csv");
    std::unique_ptr<ode::BoostSolver> solver5 = std::make_unique<ode::BoostSolver>(ode::BoostSolver(&ode, state));    
    std::vector<std::vector<double>> result5 = solver5 -> solve(0, 0.1, 100);
    save_solution(result4, file);
    file.close();
}

void ex4()
{
    std::ofstream file;
    file.open("csv/spring/euler.csv");
    ode::SpringODE ode({1, 1});
    std::vector<double> state({0, 1});
    std::unique_ptr<ode::Euler> solver = std::make_unique<ode::Euler>(ode::Euler(&ode, state));    
    std::vector<std::vector<double>> result = solver -> solve(0, 0.1, 100);
    save_solution(result, file);
    file.close();

    file.open("csv/spring/backward_euler.csv");
    std::unique_ptr<ode::BackwardEuler> solver2 = std::make_unique<ode::BackwardEuler>(ode::BackwardEuler(&ode, state));    
    std::vector<std::vector<double>> result2 = solver2 -> solve(0, 0.1, 100);
    save_solution(result2, file);
    file.close();

    file.open("csv/spring/rk2.csv");
    std::unique_ptr<ode::RungeKutta2> solver3 = std::make_unique<ode::RungeKutta2>(ode::RungeKutta2(&ode, state));    
    std::vector<std::vector<double>> result3 = solver3 -> solve(0, 0.1, 1000);
    save_solution(result3, file);
    file.close();

    file.open("csv/spring/rk4.csv");
    std::unique_ptr<ode::RungeKutta4> solver4 = std::make_unique<ode::RungeKutta4>(ode::RungeKutta4(&ode, state));    
    std::vector<std::vector<double>> result4 = solver4 -> solve(0, 0.1, 1000);
    save_solution(result3, file);
    file.close();

    file.open("csv/spring/boost.csv");
    std::unique_ptr<ode::BoostSolver> solver5 = std::make_unique<ode::BoostSolver>(ode::BoostSolver(&ode, state));    
    std::vector<std::vector<double>> result5 = solver5 -> solve(0, 0.1, 1000);
    save_solution(result4, file);
    file.close();
}

void ex5()
{
    std::ofstream file;
    file.open("csv/decay/euler.csv");
    ode::DecayODE ode({(1.0 / 25.0)});
    std::vector<double> state({1.0});
    std::unique_ptr<ode::Euler> solver = std::make_unique<ode::Euler>(ode::Euler(&ode, state));    
    std::vector<std::vector<double>> result = solver -> solve(0, 0.080001, 100);
    save_solution(result, file);
    file.close();

    file.open("csv/decay/backward_euler.csv");
    std::unique_ptr<ode::BackwardEuler> solver2 = std::make_unique<ode::BackwardEuler>(ode::BackwardEuler(&ode, state));    
    std::vector<std::vector<double>> result2 = solver2 -> solve(0, 0.080001, 100);
    save_solution(result2, file);
    file.close();

    file.open("csv/decay/rk2.csv");
    std::unique_ptr<ode::RungeKutta2> solver3 = std::make_unique<ode::RungeKutta2>(ode::RungeKutta2(&ode, state));    
    std::vector<std::vector<double>> result3 = solver3 -> solve(0, 0.080001, 100);
    save_solution(result3, file);
    file.close();

    file.open("csv/decay/rk4.csv");
    std::unique_ptr<ode::RungeKutta4> solver4 = std::make_unique<ode::RungeKutta4>(ode::RungeKutta4(&ode, state));    
    std::vector<std::vector<double>> result4 = solver4 -> solve(0, 0.080001, 100);
    save_solution(result3, file);
    file.close();

    file.open("csv/decay/boost.csv");
    std::unique_ptr<ode::BoostSolver> solver5 = std::make_unique<ode::BoostSolver>(ode::BoostSolver(&ode, state));    
    std::vector<std::vector<double>> result5 = solver5 -> solve(0, 0.080001, 100);
    save_solution(result4, file);

}

int main(void)
{
    ex2();
    ex3();
    ex4();
    ex5();
    return 0;
}
