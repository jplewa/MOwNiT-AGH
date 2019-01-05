#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include <algorithm>
#include <memory>
#include <fstream>

typedef std::vector<std::array<float, 2>> PointsList2D;

/// <summary>
/// Template function used for clamping specified datasets
/// </summary>
template <typename T>
T GetIndexClamped(const std::vector<T> points, int index)
{
	if (index < 0)
		return points[0];
	else if (index >= int(points.size()))
		return points.back();
	else
		return points[index];
}

namespace exceptions
{
	/// <summary>
	/// Helper class used for throwing not implemented error
	/// </summary>
	class NotImplementedException : public std::logic_error
	{
	public:
		NotImplementedException() : std::logic_error{ "Function not yet implemented." } {}
	};

	class IncorrectVectorsException : public std::logic_error
	{
	public:
		IncorrectVectorsException() : std::logic_error{ "Provided vectors have different lengths" } {}
	};
}

namespace interpolation
{
	/// <summary>
	/// Defines interface for interpolation classes
	/// </summary>
	class IInterpolation
	{
		// public construction and destruction methods
	public:
		virtual ~IInterpolation() = default;

		// public interface methods
	public:
		virtual void Interpolate1D(int pointsToInterpolate) = 0;
		virtual void Interpolate2D(int pointsToInterpolate) = 0;

	};

	/// <summary>
	/// Defines Hermite Cubic Interpolation 
	/// </summary>
	class CubicInterpolation : public IInterpolation
	{
		// public construction methods
	public:
		CubicInterpolation(const PointsList2D& points) : pointsList(points) {}

		// IInterpolation methods
	public:
		void Interpolate2D(int pointsToInterpolate) override
		{
			std::vector<int> index(pointsToInterpolate);
			std::vector<float> t;
			std::vector<float> tx;

			int i = 0;
			int points_size = pointsList.size() - 1;
			std::generate(index.begin(), index.end(), [&i, &pointsToInterpolate, &points_size, &t, &tx]()
			{
				float percent = ((float)i) / (float(pointsToInterpolate - 1));
				tx.push_back((points_size)* percent);	
				t.push_back(tx[i] - floor(tx[i]));
				return int(tx[i++]);
			});

			for (int i = 0; i < pointsToInterpolate; ++i)
			{
				std::array<PolynomialCoeffs, 2> coeffs;
				std::array<float, 2> A = GetIndexClamped(pointsList, index[i] - 1);
				std::array<float, 2> B = GetIndexClamped(pointsList, index[i]);
				std::array<float, 2> C = GetIndexClamped(pointsList, index[i] + 1);
				std::array<float, 2> D = GetIndexClamped(pointsList, index[i] + 2);

				for (int j = 0; j < 2; j++)
				{
					coeffs[j].A = A[j];
					coeffs[j].B = B[j];
					coeffs[j].C = C[j];
					coeffs[j].D = D[j];
				}

				float x = CubicHermite(coeffs[0], t[i]);
				float y = CubicHermite(coeffs[1], t[i]);

				std::cout << "Value at " << tx[i] << " = " << x << "  " << y << std::endl;
			}
		};

		void Interpolate1D(int pointsToInterpolate) override
		{
			std::ofstream outfile;
			outfile.open("csv/cubic.csv", std::ios::out);

			std::vector<int> index(pointsToInterpolate);
			std::vector<float> t;
			std::vector<float> tx;

			int i = 0;
			int points_size = pointsList.size() - 1;
			std::generate(index.begin(), index.end(), [&i, &pointsToInterpolate, &points_size, &t, &tx]()
			{
				float percent = ((float) i) / (float(pointsToInterpolate - 1));
				tx.push_back((points_size)* percent);	
				t.push_back(tx[i] - floor(tx[i]));
				return int(tx[i++]);
			});

			for (int i = 0; i < pointsToInterpolate; ++i)
			{
				PolynomialCoeffs coeffs;
				std::array<float, 2> A = GetIndexClamped(pointsList, index[i] - 1);
				std::array<float, 2> B = GetIndexClamped(pointsList, index[i]);
				std::array<float, 2> C = GetIndexClamped(pointsList, index[i] + 1);
				std::array<float, 2> D = GetIndexClamped(pointsList, index[i] + 2);

				coeffs.A = A[0];
				coeffs.B = B[0];
				coeffs.C = C[0];
				coeffs.D = D[0];
				
				float x = CubicHermite(coeffs, t[i]);	

				std::cout << "Value at " << tx[i] << " = " << x << std::endl;
				outfile << tx[i] << "," << x << std::endl;
			}
			outfile.close();
		};

	// private methods
	private:
		struct PolynomialCoeffs
		{
			float A, B, C, D, t;
		};

		float CubicHermite(PolynomialCoeffs coeffs, float t) const
		{
			float a = -coeffs.A / 2.0f + (3.0f*coeffs.B) / 2.0f - (3.0f*coeffs.C) / 2.0f + coeffs.D / 2.0f;
			float b = coeffs.A - (5.0f*coeffs.B) / 2.0f + 2.0f*coeffs.C - coeffs.D / 2.0f;
			float c = -coeffs.A / 2.0f + coeffs.C / 2.0f;
			float d = coeffs.B;

			return a * pow(t, 3) + b * pow(t, 2) + c * t + d;
		}

	// private members
	private:
		const PointsList2D& pointsList;
	};

	class LagrangeInterpolation : public IInterpolation
	{
	public:
		LagrangeInterpolation(const PointsList2D& points) : pointsList(points) {}

		void Interpolate1D(int pointsToInterpolate) override
		{			
			std::ofstream outfile;
			outfile.open("csv/lagrange.csv", std::ios::out);
			for (int i = 0; i < pointsToInterpolate; ++i)
			{
				float percent = ((float)i) / (float(pointsToInterpolate - 1));
            	float x = pointsList.back()[0] * percent;
            	float y = LagrangeInterpolate(x);
				
				std::cout << "Value at " << x << " = " << y << std::endl;
				outfile << x << "," << y << std::endl;
			}
			outfile.close();
		};

		void Interpolate2D(int pointsToInterpolate) override
		{
		};

	private:

		float LagrangeInterpolate (float x)
		{
		    float sum = 0.0f;
		    for (int i = 0; i < pointsList.size(); ++i)
			{
				float lx = 1.0;
				for (int j = 0; j < pointsList.size(); ++j) 
				{
    		    	if (j != i)
					{
    		        	lx *= (x - pointsList[j][0]) / (pointsList[i][0] - pointsList[j][0]);
					}
    			}
		        sum += pointsList[i][1] * lx;
			}
		    return sum;
		}

	private:
		const PointsList2D& pointsList;
	};
}

namespace regression
{
	float mse(std::vector<float> x, std::vector<float> y)
	{
		if (x.size() != y.size()){
			throw exceptions::IncorrectVectorsException();
		}
		float mse = 0;
		for (int i = 0; i < x.size(); i++)
		{
			mse += pow(x[i] - y[i], 2);
		}
		return mse / x.size();
	}

	float vectorElemSum(std::vector<float> v)
	{
		float result = 0;
		for (int i = 0; i < v.size(); i++)
		{
			result += v[i];
		}
		return result;
	}

	float vectorElemMean(std::vector<float> v)
	{
		return vectorElemSum(v) / v.size();
	}

	std::vector<float> subtractValueFromVectorElems(std::vector<float> x, float y)
	{
		std::vector<float> result;
		for (int i = 0; i < x.size(); i++)
		{
			result.push_back(x[i] - y);
		}
		return result;
	}	
	std::vector<float> multiplyVectorElems(std::vector<float> x, std::vector<float> y)
	{
		std::vector<float> result;
		for (int i = 0; i < x.size(); i++)
		{
			result.push_back(x[i] * y[i]);
		}
		return result;
	}

	std::pair<float, float> linearRegression(std::vector<float> x, std::vector<float> y)
	{
		float xMean = vectorElemMean(x);
		float yMean = vectorElemMean(y);

		std::vector<float> xMinusMean = subtractValueFromVectorElems(x, xMean);
		std::vector<float> yMinusMean = subtractValueFromVectorElems(y, yMean);

		float a = vectorElemSum(multiplyVectorElems(xMinusMean, yMinusMean));
		a /= vectorElemSum(multiplyVectorElems(xMinusMean, xMinusMean));
		float b = yMean - a * xMean;

		std::cout << "y = " << a << " * x + "  << b << std::endl;

		return std::make_pair(a, b);
	}

	class LinearRegressor
	{
	    public:
	        
			LinearRegressor() { }
	        
			float fit(std::vector<float> x, std::vector<float> y) 
	        { 
				coefficients = regression::linearRegression(x, y);
				float mse = regression::mse(y, predict(x));

	            return mse;
	        }
	        
			std::vector<float> predict(std::vector<float> x) 
	        { 
				std::vector<float> y;
				for (int i = 0; i < x.size(); i++)
				{
					y.push_back(coefficients.first * x[i] + coefficients.second);
				}
	            return y; 
	        }
	
	    private:
	        std::pair<float, float> coefficients;
	};
}

void ex1()
{
	const PointsList2D ex1Points =
	{
		{ 0.0f, 1.1f },
		{ 1.6f, 8.3f },
		{ 2.3f, 6.5f },
		{ 3.5f, 4.7f },
		{ 4.3f, 3.1f },
		{ 5.9f, 7.5f },
		{ 6.8f, 0.0f },
	};
	std::unique_ptr<interpolation::CubicInterpolation> cubicInterpolation = 
					std::make_unique<interpolation::CubicInterpolation>(interpolation::CubicInterpolation(ex1Points));
	std::cout << "1D cubic interpolation" << std::endl;
	cubicInterpolation -> Interpolate1D(13);
	std::cout << std::endl;
}

void ex2()
{
	const PointsList2D ex2Points =
	{
		{ 0.0f, 1.1f },
		{ 1.6f, 8.3f },
		{ 2.3f, 6.5f },
		{ 3.5f, 4.7f },
		{ 4.3f, 3.1f },
		{ 5.9f, 7.5f },
		{ 6.8f, 0.0f },
	};
	std::unique_ptr<interpolation::LagrangeInterpolation> lagrangeInterpolation = 
					std::make_unique<interpolation::LagrangeInterpolation>(interpolation::LagrangeInterpolation(ex2Points));
	std::cout << "1D Lagrange interpolation" << std::endl;
	lagrangeInterpolation -> Interpolate1D(13);
	std::cout << std::endl;
}

void ex2Comparison()
{
	const PointsList2D ex2ComparisonCubic =
	{
		{ 0.0f, 1.1f }, { 1.6f, 8.3f }, { 2.3f, 6.5f }, { 3.5f, 4.7f }, { 4.3f, 3.1f },
		{ 5.9f, 7.5f }, { 6.8f, 0.0f }, { 1.2f, 7.2f }, { 8.1f, 8.1f }, { 20.5f, 20.5f },
		{ 10.0f, 10.0f }, { 1.8f, 1.8f }, { 5.2f, 5.2f }, { 10.1f, 10.1f }, { 3.5f, 3.5f },
		{ 8.0f, 15.0f }, { 1.8f, 16.0f }, { 5.2f, 17.0f },{ 10.1f, 18.0f },	{ 3.5f, 19.0f },
		{ 8.0f, 20.0f },
	};
	std::unique_ptr<interpolation::CubicInterpolation> cubicInterpolation = 
					std::make_unique<interpolation::CubicInterpolation>(interpolation::CubicInterpolation(ex2ComparisonCubic));
	std::cout << "1D cubic interpolation" << std::endl;
	cubicInterpolation -> Interpolate1D(10000);
	std::cout << std::endl;

	const PointsList2D ex2ComparisonLagrange =
	{
		{ 0.0f,  0.0f }, { 1.0f,  1.6f }, { 2.0f,  2.3f }, { 3.0f,  3.5f },  { 4.0f,  4.3f },
		{ 5.0f,  5.9f }, { 6.0f,  6.8f }, { 7.0f,  1.2f },  { 8.0f,  8.1f }, { 9.0f,  20.5f },
		{ 10.0f, 10.0f }, { 11.0f,  1.8f }, { 12.0f,  5.2f }, { 13.0f,  10.1f }, { 14.0f,  3.5f },
		{ 15.0f,  8.0f }, { 16.0f,  1.8f }, { 17.0f,  5.2f }, { 18.0f,  10.1f }, { 19.0f,  3.5f },
		{ 20.0f,  8.0f },
	};
	std::unique_ptr<interpolation::LagrangeInterpolation> lagrangeInterpolation = 
					std::make_unique<interpolation::LagrangeInterpolation>(interpolation::LagrangeInterpolation(ex2ComparisonLagrange));
	std::cout << "1D Lagrange interpolation" << std::endl;
	lagrangeInterpolation -> Interpolate1D(10000);
	std::cout << std::endl;
}

void ex5()
{
	std::vector<float> xs {108, 19, 13, 124, 40, 57, 23, 14, 45, 10,
						   5, 48, 11, 23, 7, 2, 24, 6, 3, 23, 6, 9,
						   9, 3, 29, 7, 4, 20,7, 4, 0, 25, 6, 5, 22, 
						   11, 61, 12, 4, 16, 13, 60, 41, 37, 55, 41,
						   11, 27, 8, 3, 17, 13, 13, 15, 8, 29, 30,
						   24, 9, 31, 14, 53, 26};
	std::vector<float> ys {392.5, 46.2, 15.7, 422.2, 119.4, 170.9, 56.9, 
						   77.5, 214, 65.3, 20.9, 248.1, 23.5, 39.6, 48.8, 
						   6.6, 134.9, 50.9, 4.4, 113, 14.8, 48.7, 52.1, 
						   13.2, 103.9, 77.5, 11.8, 98.1, 27.9, 38.1, 
						   0, 69.2, 14.6, 40.3, 161.5, 57.2, 217.6, 
						   58.1, 12.6, 59.6, 89.9, 202.4, 181.3, 
						   152.8, 162.8, 73.4, 21.3, 92.6, 76.1, 
						   39.9, 142.1, 93, 31.9, 32.1, 55.6, 133.3, 
						   194.5, 137.9, 87.4, 209.8, 95.5, 244.6, 187.5};

	std::cout << "Linear regression" << std::endl;
	std::pair<float, float> ab = regression::linearRegression(xs, ys);
	std::cout << std::endl;
}

void ex6()
{
	std::vector<float> xs {108, 19, 13, 124, 40, 57, 23, 14, 45, 10,
						   5, 48, 11, 23, 7, 2, 24, 6, 3, 23, 6, 9,
						   9, 3, 29, 7, 4, 20,7, 4, 0, 25, 6, 5, 22, 
						   11, 61, 12, 4, 16, 13, 60, 41, 37, 55, 41,
						   11, 27, 8, 3, 17, 13, 13, 15, 8, 29, 30,
						   24, 9, 31, 14, 53, 26};
	std::vector<float> ys {392.5, 46.2, 15.7, 422.2, 119.4, 170.9, 56.9, 
						   77.5, 214, 65.3, 20.9, 248.1, 23.5, 39.6, 48.8, 
						   6.6, 134.9, 50.9, 4.4, 113, 14.8, 48.7, 52.1, 
						   13.2, 103.9, 77.5, 11.8, 98.1, 27.9, 38.1, 
						   0, 69.2, 14.6, 40.3, 161.5, 57.2, 217.6, 
						   58.1, 12.6, 59.6, 89.9, 202.4, 181.3, 
						   152.8, 162.8, 73.4, 21.3, 92.6, 76.1, 
						   39.9, 142.1, 93, 31.9, 32.1, 55.6, 133.3, 
						   194.5, 137.9, 87.4, 209.8, 95.5, 244.6, 187.5};

	std::cout << "Linear regression" << std::endl;

	std::unique_ptr<regression::LinearRegressor> linearRegression = 
					std::make_unique<regression::LinearRegressor>(regression::LinearRegressor());

	float mse = linearRegression -> fit(xs, ys);

	std::cout << "MSE(y, predicted_y) = " << mse << std::endl;
	std::cout << "RMSE(y, predicted_y) = " << sqrt(mse) << std::endl;

	std::vector<float> predictedPoints = linearRegression -> predict(xs);

	for(int i = 0; i < predictedPoints.size(); i++){
		std::cout << "f(" << xs[i] << ") = " << ys[i] << " -> " << predictedPoints[i] << std::endl;
	}
}

int main()
{
	//ex1();
	//ex2();
	ex2Comparison();
	//ex5();
	//ex6();
	
	return 0;
}