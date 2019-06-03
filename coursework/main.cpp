#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <thread>
#include <future>
#include "lib.h"
#include "fem.h"

using namespace std;
using namespace placeholders;

void investigate_t_changing(
	int n,
	const string& filename,
	const function<double(double)>& ft
) {
	double default_value = ft(0);

	grid_generator_t grid(-1, 1, n);
	vector<future<pair<double, double>>> results;
	for (int i = 0; i < grid.size(); i++) {
		results.push_back(async([i, ft, grid] () -> pair<double, double> {
			double res;
			double time = calc_time_microseconds([&] () {
				res = ft(grid(i));
			});
			return {res, time};
		}));
	}

	int counter = 0;
	double max = -numeric_limits<double>::max();
	vector<pair<double, double>> matrix;
	for (auto& i : results) {
		write_percent(double(counter)/results.size());
		auto value = i.get();
		matrix.push_back(value);
		if (value.first == value.first)
			if (value.first > max)
				max = value.first;
		counter++;
	}
	cout << "\r                \r";

	ofstream fout(filename + ".txt");
	fout << "t\tres\ttime" << endl;
	for (auto& i : matrix) {
		auto value = i;
		if (value.first != value.first) value.first = max;
		fout << grid(counter) << "\t" << value.first << "\t" << value.second << endl;
	}
	fout.close();
}

void investigate_t2_changing(
	int n,
	const string& filename,
	const function<double(double, double)>& ft
) {
	double default_value = ft(0, 0);

	grid_generator_t grid(-1, 1, n);
	vector<vector<future<pair<double, double>>>> results;
	for (int i = 0; i < grid.size(); i++) {
		results.push_back({});
		for (int j = 0; j < grid.size(); j++) {
			results.back().push_back(async([i, j, ft, grid] () -> pair<double, double> {
				double res;
				double time = calc_time_microseconds([&] () {
					res = ft(grid(i), grid(j));
				});
				return {res, time};
			}));
		}
	}

	int counter = 0;
	double max = -numeric_limits<double>::max();
	vector<vector<pair<double, double>>> matrix;
	for (auto& i : results) {
		write_percent(double(counter)/results.size());
		matrix.push_back({});
		for (auto& j : i) {
			auto value = j.get();
			matrix.back().push_back(value);
			if (value.first == value.first)
				if (value.first > max)
					max = value.first;
		}
		counter++;
	}
	cout << "\r                \r";

	ofstream fout(filename + ".data.txt");
	ofstream fout2(filename + ".time.txt");
	for (auto& i : matrix) {
		for (auto& j : i) {
			auto value = j;
			if (value.first != value.first) value.first = max;
			fout << value.first << "\t";
			fout2 << value.second << "\t";
		}
		fout << endl;
		fout2 << endl;
	}
	fout.close();
	fout2.close();

	fout.open(filename + ".x.txt");
	for (int i = 0; i < grid.size(); i++)
		fout << grid(i) << endl;
	fout.close();

	fout.open(filename + ".y.txt");
	for (int i = 0; i < grid.size(); i++)
		fout << grid(i) << endl;
	fout.close();
}

int main() {
	/*constants_t c = {1, 1, 1, 1};
	//auto u = [] (double x, double y, double t) -> double { return exp(x*x*(y-t)) + t*t; };
	auto u = [] (double x, double y, double t) -> double { return x*x + y*y + t*t; };
	auto f = calc_right_part(u, c);
	boundary_setter_t set_boundary_conditions = bind(write_first_boundary_conditions, _1, _2, _3, _4, u);

	grid_t grid;
	grid.calc(grid_generator_t(0, 1, 1), grid_generator_t(0, 1, 1));

	grid_generator_t time_grid(0, 1, 1, 0);

	vector_t q0 = calc_true_approx(bind(u, _1, _2, time_grid(0)), grid.bes);
	vector_t q1 = calc_true_approx(bind(u, _1, _2, time_grid(1)), grid.bes);
	vector_t q = calc_true_approx(bind(u, _1, _2, time_grid(time_grid.size()-1)), grid.bes);

	auto res = solve_differential_equation(f, set_boundary_conditions, q0, q1, c, grid, time_grid);

	double norm = (q-res.back()).norm();
	cout << "norm: " << norm << endl;

	double norm_integral = calc_integral_norm(bind(u, _1, _2, time_grid(time_grid.size() - 1)), grid.es, res.back());
	cout << "norm integral: " << norm_integral << endl;

	system("pause");*/
	cout << calc_time_microseconds([](){
		investigate_t2_changing(
			100,
			"space_nonlinear_grid",
			[] (double tx, double ty) -> double {
				constants_t c = {1, 1, 1, 1};
				auto u = [] (double x, double y, double t) -> double { return exp(x*y) + exp(t*t); };
				auto f = calc_right_part(u, c);
				boundary_setter_t set_boundary_conditions = bind(write_first_boundary_conditions, _1, _2, _3, _4, u);

				grid_t grid;
				grid.calc(grid_generator_t(0, 1, 10, tx), grid_generator_t(0, 1, 10, ty));

				grid_generator_t time_grid(0, 1, 10);

				vector_t q0 = calc_true_approx(bind(u, _1, _2, time_grid(0)), grid.bes);
				vector_t q1 = calc_true_approx(bind(u, _1, _2, time_grid(1)), grid.bes);

				auto res = solve_differential_equation(f, set_boundary_conditions, q0, q1, c, grid, time_grid);

				double norm_integral = calc_integral_norm(bind(u, _1, _2, time_grid(time_grid.size() - 1)), grid.es, res.back());
				return norm_integral;
			}
		);

		investigate_t_changing(
			1000,
			"time_nonlinear_grid",
			[] (double t) -> double {
				constants_t c = {1, 1, 1, 1};
				auto u = [] (double x, double y, double t) -> double { return exp(x*y) + exp(t*t); };
				auto f = calc_right_part(u, c);
				boundary_setter_t set_boundary_conditions = bind(write_first_boundary_conditions, _1, _2, _3, _4, u);

				grid_t grid;
				grid.calc(grid_generator_t(0, 1, 10), grid_generator_t(0, 1, 10));

				grid_generator_t time_grid(0, 1, 10, t);

				vector_t q0 = calc_true_approx(bind(u, _1, _2, time_grid(0)), grid.bes);
				vector_t q1 = calc_true_approx(bind(u, _1, _2, time_grid(1)), grid.bes);

				auto res = solve_differential_equation(f, set_boundary_conditions, q0, q1, c, grid, time_grid);

				double norm_integral = calc_integral_norm(bind(u, _1, _2, time_grid(time_grid.size() - 1)), grid.es, res.back());
				return norm_integral;
			}
		);
	})/1000/1000 << "s" << endl;
	system("pause");
}