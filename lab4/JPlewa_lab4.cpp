#include <algorithm>
#include <array>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <numeric>
#include <vector>
#include <fftw3.h>
#include <iostream>
#include <chrono>
#include <limits>
#include <string>

#define CUT_OFF_SIZE 2

using namespace std::complex_literals;  

namespace exceptions
{
    class IncorrectVectorSizeException : public std::logic_error
    {
    public:
        IncorrectVectorSizeException() : std::logic_error{"Vector size must be a power of 2"} {}
    };
}

namespace fourier
{
    template<typename T>
    class ITransform
    {
    public:
        virtual ~ITransform() = default;

    public:
        virtual std::vector <std::complex<T>> calculate() = 0;

    };

    template<typename T>
    class DFT : public ITransform<T>
    {

    public:
        DFT(std::vector <std::complex<T>> &input) : rawDataset(input) {}

    public:
        std::vector <std::complex<T>> calculate() override
        {
            // prepare result vector
            std::vector <std::complex<T>> transformedDataset(rawDataset.size(), std::complex<T>(0, 0));
            // compute result vector elements one by one
            for (int k = 0; k < rawDataset.size(); k++)
            {
                // each result vector element is a sum
                for (int n = 0; n < rawDataset.size(); n++)
                {
                    // using Euler's formula
                    std::complex <T> tmp(cos(2 * M_PI * k * n / (T) rawDataset.size()), -(sin(2 * M_PI * k * n / (T) rawDataset.size())));
                    transformedDataset[k] += rawDataset[n] * tmp;
                }
            }
            return transformedDataset;
        }

    private:
        const std::vector <std::complex<T>> &rawDataset;
    };

    
    template<typename T>
    class FFT : public ITransform<T>
    {

    public:
        FFT(std::vector <std::complex<T>> input) : rawDataset(input) {}


    private:
        // recursive inner function
        void cooleyTurkey(std::vector <std::complex<T>> &transformedDataset, int index, int begin, int N, int step)
        {
            // base case DFT
            if (N <= CUT_OFF_SIZE)
            {
                // copy data into a new vector
                std::vector <std::complex<T>> sliceCopy(N);
                for (int i = 0; i < N; i++)
                {
                    sliceCopy[i] = rawDataset[begin + i * step];
                }
                // calculate the DFT
                std::unique_ptr <fourier::DFT<double>> dft = std::make_unique<fourier::DFT<double>>(fourier::DFT<double>(sliceCopy));
                std::vector <std::complex<double>> dft_result = dft -> calculate();
                // copy results back into the result vector
                for (int i = 0; i < N; i++)
                {
                    transformedDataset[index + i] = dft_result[i];
                }
            }
            else
            {
                // calculate left- and right-half FFT
                // each sub-vector has half less elements and includes every other element from the currect vector 
                // (with an offset in case of the right-half vector)
                cooleyTurkey(transformedDataset, index,         begin,        N / 2, 2 * step);
                cooleyTurkey(transformedDataset, index + N / 2, begin + step, N / 2, 2 * step);
                // use symmetry to calculate final results
                for (int k = index; k < index + (N / 2); k++)
                {
                    std::complex<T> t = transformedDataset[k];
                    std::complex<T> tmp(cos(2 * M_PI * k / N), -sin(2 * M_PI * k / N));
                    transformedDataset[k]         = t + transformedDataset[k + N / 2] * tmp;
                    transformedDataset[k + N / 2] = t - transformedDataset[k + N / 2] * tmp;
                }
            }
        }

    public:

        std::vector <std::complex<T>> calculate() override
        {
            if (ceil(log2(rawDataset.size())) != floor(log2(rawDataset.size())))
            {
                throw exceptions::IncorrectVectorSizeException();
            }
            std::vector <std::complex<T>> transformedDataset(rawDataset.size());
            // out-of-place radix-2 DIT (decimation-in-time) FFT
            cooleyTurkey(transformedDataset, 0, 0, rawDataset.size(), 1);
            return transformedDataset;
        }

    private:
        const std::vector <std::complex<T>> rawDataset;
    };
}

template<typename T>
std::ostream &operator<<(std::ostream &output, std::vector <T> const &values)
{
    for (auto const &value : values)
    {
        output << value << " ";
    }
    output << std::endl;
    return output;
}

void initializeVector(std::vector <std::complex<double>> &fVector)
{
    srand(time(nullptr));
    for (int i = 0; i < fVector.size(); i++)
    {
        fVector[i] = rand() / double(RAND_MAX) + i * rand() / double(RAND_MAX);
    }
}

void measureDFT(int benchmarkLoops, std::vector <std::complex<double>> fVector, std::ofstream &times)
{
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    std::unique_ptr <fourier::DFT<double>> dft = std::make_unique<fourier::DFT<double>>(fourier::DFT<double>(fVector));
    
    long result = 0;
    std::vector <std::complex<double>> dft_result;

    for (int i = 0; i < benchmarkLoops; i ++)
    {
        begin = std::chrono::steady_clock::now();
        dft_result = dft -> calculate();
        end = std::chrono::steady_clock::now();
        result += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    }
    times << ((double) result / benchmarkLoops) / 1000;
}

void measureFFT(int benchmarkLoops, std::vector <std::complex<double>> fVector, std::ofstream &times)
{
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    std::unique_ptr <fourier::FFT<double>> fft = std::make_unique<fourier::FFT<double>>(fourier::FFT<double>(fVector));
    
    long result = 0;
    std::vector <std::complex<double>> fft_result;

    for (int i = 0; i < benchmarkLoops; i ++)
    {
        begin = std::chrono::steady_clock::now();
        fft_result = fft -> calculate();
        end = std::chrono::steady_clock::now();
        result += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    }
    times << ((double) result / benchmarkLoops) / 1000 << ",";
}

void ex1()
{
    std::vector <std::complex<double>> dft_example_1{1, 2-1i, -1i, -1+2i};
    std::unique_ptr <fourier::DFT<double>> dft = std::make_unique<fourier::DFT<double>>(fourier::DFT<double>(dft_example_1));
    std::vector <std::complex<double>> dft_result = dft -> calculate();
    std::cout << dft_result << std::endl;
}

void ex2()
{
    std::vector <std::complex<double>> fft_example_1{1, 2-1i, -1i, -1+2i};
    std::unique_ptr <fourier::FFT<double>> fft = std::make_unique<fourier::FFT<double>>(fourier::FFT<double>(fft_example_1));
    std::vector <std::complex<double>> fft_result = fft -> calculate();
    std::cout << fft_result << std::endl;
}

void ex3()
{
    int N = 1;
    std::ofstream times;
    times.open("csv/ex3.csv");
    for (int i = 0; i < 12; i++)
    {   
        for (int j = 0; j < 10; j++)
        {
            times << N << ",";
            std::vector <std::complex<double>> fVector(N);
            initializeVector(fVector);
            measureFFT(10, fVector, times);
            measureDFT(10, fVector, times);
            times << std::endl;
        }
        N *= 2;
    }
    times.close();

}

void ex4()
{
    std::vector <std::complex<double>> shampoo;

    std::string line;
    std::ifstream in("csv/ex4_in.csv");
    if (in.is_open())
    {
        while (getline(in, line))
        {
            std::string month;
            double value;
            getline(in, month, ',');
            in >> value;
            shampoo.push_back(std::complex<double>(value, 0));
        }
        in.close();
    }
    else std::cout << "Unable to open file"; 
    std::unique_ptr <fourier::FFT<double>> fft = std::make_unique<fourier::FFT<double>>(fourier::FFT<double>(shampoo));
    std::vector <std::complex<double>> fft_result = fft -> calculate();
    std::ofstream out("csv/ex4_out.csv");
    if (out.is_open())
    {
        for (int i = 0; i < fft_result.size(); i++)
        {
            out << fft_result[i].real() << "," << fft_result[i].imag() << std::endl;
        }
        out.close();
    }
}

void measureFFTW(int benchmarkLoops, std::vector <std::complex<double>> fVector, std::ofstream &times)
{
    int N = fVector.size();
    fftw_complex* in = new fftw_complex[N];
    fftw_complex* out = new fftw_complex[N];
    
    fftw_plan p, q;
    for (int i = 0; i < N; i++) {
        in[i][0] = fVector[i].real();
        in[i][1] = fVector[i].imag();
    }
    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;
   
    long result = 0;

    for (int i = 0; i < benchmarkLoops; i ++)
    {
        begin = std::chrono::steady_clock::now();
        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);       
         end = std::chrono::steady_clock::now();
        result += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
    }
    times << ((double) result / benchmarkLoops) / 1000 << ",";
}


void ex5()
{
    int N = 1;
    std::ofstream times;
    times.open ("csv/ex5.csv");
    
    for (int i = 0; i < 20; i++)
    {   
        times << N << ",";
        std::vector <std::complex<double>> fVector(N);
        initializeVector(fVector);
        measureFFT(10, fVector, times);
        measureFFTW(10, fVector, times);
        times << std::endl;
        N *= 2;
    }
}

int main()
{
    ex1();
    ex2();
    ex3();
    ex4();
    ex5();
    return 0;
}
