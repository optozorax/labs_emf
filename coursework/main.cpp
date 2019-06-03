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

template<class Ret, class Key>
class async_performer_t
{
public:
	void add(const function<Ret(void)>& f, const Key& key) {
		mf[key] = async(f);
	}

	void finish(void) {
		int counter = 0;
		for (auto& i : mf) {
			if (counter % (mf.size()/10000 + 1) == 0)
				write_percent(double(counter)/mf.size());
			auto value = i.second.get();
			m[i.first] = value;
			counter++;
		}
		cout << "\r       \r";
	}

	auto begin(void) { return m.begin(); }
	auto end(void) { return m.end(); }
private:
	map<Key, future<Ret>> mf;
	map<Key, Ret> m;
};

void investigate_t_changing(
	int n,
	const string& filename,
	const function<double(double)>& ft
) {
	double default_value = ft(0);

	async_performer_t<pair<double, double>, int> performer;

	grid_generator_t grid(-1, 1, n);
	for (int i = 0; i < grid.size(); i++) {
		performer.add([i, ft, grid] () -> pair<double, double> {
			double res;
			double time = calc_time_microseconds([&] () {
				res = ft(grid(i));
			});
			return {res, time};
		}, i);
	}

	performer.finish();

	int counter = 0;
	auto mymax = max_element(performer.begin(), performer.end(), [] (auto& a, auto& b) -> bool { 
		if (isnan(a.second.first)) {
			return true;
		} else {
			return a.second.first < b.second.first; 
		}
	});

	ofstream fout(filename + ".txt");
	fout << "t\tres\ttime" << endl;
	for (auto& i : performer) {
		fout 
			<< grid(i.first) << "\t" << 
			(isnan(i.second.first) ? mymax->second.first : i.second.first) << "\t" << 
			i.second.second << endl;
	}
	fout.close();
}

void investigate_t2_changing(
	int n,
	const string& filename,
	const function<double(double, double)>& ft
) {
	double default_value = ft(0, 0);

	async_performer_t<pair<double, double>, pair<int, int>> performer;

	grid_generator_t grid(-1, 1, n);
	for (int i = 0; i < grid.size(); i++) {
		for (int j = 0; j < grid.size(); j++) {
			performer.add([i, j, ft, grid] () -> pair<double, double> {
				double res;
				double time = calc_time_microseconds([&] () {
					res = ft(grid(i), grid(j));
				});
				return {res, time};
			}, {i, j});
		}
	}

	performer.finish();

	int counter = 0;
	auto mymax = max_element(performer.begin(), performer.end(), [] (auto& a, auto& b) -> bool { 
		if (isnan(a.second.first)) {
			return true;
		} else {
			return a.second.first < b.second.first; 
		}
	});

	ofstream fout(filename + ".data.txt");
	ofstream fout2(filename + ".time.txt");
	int last_line = 0;
	for (auto& i : performer) {
		if (last_line != i.first.first) {
			fout << endl;
			fout2 << endl;
		}
		fout << (isnan(i.second.first) ? mymax->second.first : i.second.first) << "\t";
		fout2 << i.second.second << "\t";
		last_line = i.first.first;
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
		vector<pair<function_3d_t, int>> u_space_mas;
		u_space_mas.push_back({[] (double x, double y, double t) -> double { return exp(x*y) + exp(t*t); }, 0});
		u_space_mas.push_back({[] (double x, double y, double t) -> double { return exp((1-x)*(1-y)) + exp((1-t)*(1-t)); }, 1});
		u_space_mas.push_back({[] (double x, double y, double t) -> double { return x*x*x + y*y*y*y*x*x*t + t*t*exp(t); }, 2});

		for (auto& i : u_space_mas) {
			auto& u = i.first;
			investigate_t2_changing(
				100,
				"space_tgrid_" + to_string(i.second),
				[u] (double tx, double ty) -> double {
					constants_t c = {1, 1, 1, 1};
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
				"time_tgrid_" + to_string(i.second),
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
		}
	})/1000/1000 << "s" << endl;
	system("pause");
}