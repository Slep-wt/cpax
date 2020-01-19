#define _USE_MATH_DEFINES

#include <jsoncpp.cpp>
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <string>
#include <fstream>
#include <cmath>

// Utilities and other stuff
std::ofstream fs;

class utils {
public:
	static std::string createOutput(const std::vector<double> &input) {
		std::string out = "[";
		for (unsigned int i = 0; i < (input.size() - 1); i++) {
			out += std::to_string(input[i]) + ", ";
		}
		out += std::to_string(input.back()) + "]";
		return out;
	}

	static void writeToFile(const std::vector<double> &alpha, const std::vector<double> &num, const int &lit) {
		fs.open("values_l" + std::to_string(lit) + ".dat", std::ios::trunc);
		for (unsigned int i = 0; i < alpha.size(); i++) {
			fs << alpha[i] << " " << num[i] << std::endl;
		}
		fs.close();
	}
};

static int run()
{
	Json::Value jres;
	std::ifstream j("settings.json", std::ifstream::binary);
	j >> jres;
	// PARAMS

	const double vc = jres["vc"].asDouble(); //2
	const double dt = jres["dt"].asDouble(); //0.1
	const double alpha = jres["alpha"].asDouble(); //0.5 
	const double alpha2inc = jres["alpha2inc"].asDouble();//0.001
	const double alpha2min = jres["alpha2min"].asDouble(); //0.1

	const double phi0 = 2 * M_PI;
	const int n = jres["n"].asDouble(); //100
	const int N = jres["N"].asDouble(); //10000

	double sens = jres["sens"].asDouble(); //0.0

	double level = jres["level"].asDouble();//2.0
	const double lin = jres["linc"].asDouble(); //2.0
	const int lit = jres["lit"].asDouble(); //1
	const double frad = jres["rad"].asDouble();
	//const double pos[2] = { jres["pos"]["x"].asDouble() , jres["pos"]["y"].asDouble() };

	// END PARAMS

	std::vector<double> minalphas;
	std::vector<double> minpoints;
	std::default_random_engine generator;
	generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	auto et_start = std::chrono::high_resolution_clock::now();
	for (unsigned int k = 1; k < lit + 1; k++)
	{
		double alpha2 = alpha2min;
		std::vector<double> steparr;
		std::vector<double> alphas;

		double rad = sqrt(log(level) / alpha);
		std::cout << "Level: " << level << std::endl;
		std::vector<double> points;

		for (unsigned int j = 0; j < n; j++) {
			std::cout << "Alpha2: " << alpha2 << std::endl;
			double savg = 0.0;
#pragma omp parallel for
			for (int i = 0; i < N; i++) {
				std::vector<double> x;
				std::vector<double> y;
				double rpa = phi0 * ((double)rand() / RAND_MAX);
				double pos[3] = { frad * cos(rpa), frad * sin(rpa), frad * sin(rpa) };

				x.push_back(pos[0]);
				y.push_back(pos[1]);

				double philast = 0.0;
				double phi = phi0 * ((double)rand() / RAND_MAX);
				x.push_back(x[0] + vc * cos(phi)*dt);
				y.push_back(y[0] + vc * sin(phi)*dt);
				unsigned int steps = 1;
				double prevconc = exp(-alpha * (pow(x[0], 2.0) + pow(y[0], 2.0)));
				while (sqrt(pow(x[steps], 2.0) + pow(y[steps], 2.0)) > rad) {
					double conc = exp(-alpha * (pow(x[steps], 2.0) + pow(y[steps], 2.0)));
					if (sens >= abs(conc - prevconc)) {
						phi = phi0 * ((double)rand() / RAND_MAX);
					}
					else if (conc > prevconc) {
						std::normal_distribution<> distribution(philast, alpha2);
						phi = distribution(generator);
					}
					else {
						std::normal_distribution<> distribution(philast + M_PI, alpha2);
						phi = distribution(generator);
					}
					x.push_back(x[steps] + vc * cos(phi)*dt);
					y.push_back(y[steps] + vc * sin(phi)*dt);
					philast = phi;
					prevconc = conc;
					steps++;
				}
				savg += steps;
			}
			steparr.push_back(savg / N);

			alphas.push_back(alpha2);

			alpha2 += alpha2inc;
		}
		std::cout << "Values (Alphas): " << utils::createOutput(alphas) << std::endl;
		std::cout << "Values (Steps) : " << utils::createOutput(steparr) << "\n" << std::endl;

		utils::writeToFile(alphas, steparr, k);
		minalphas.push_back(alphas.at(distance(begin(steparr), min_element(begin(steparr), end(steparr)))));
		level = pow(level, lin);
	}

	fs.open("minalphas.dat", std::ios::trunc);

	for (unsigned int i = 0; i < minalphas.size(); i++) {
		fs << minalphas.at(i) << std::endl;
	}
	auto et_stop = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - et_start).count();
	std::cout << "Execution time: " << et_stop << " ms" << "\n" << "Press ANY key to exit...";

	std::cin.get();
	return 0;
}

int _stdcall main() {
	return run();
}
