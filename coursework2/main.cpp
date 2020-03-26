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

//-----------------------------------------------------------------------------
struct fem_result_t
{
	double integral_residual;
	double norm_residual;
	double time;
};

//-----------------------------------------------------------------------------
fem_result_t calc_fem_residual(
	const function_3d_t& u,
	const grid_generator_t& x_grid,
	const grid_generator_t& y_grid,
	const grid_generator_t& time_grid,
	const constants_t& c = {1, 1, 1, 1}
) {
	fem_result_t res;
	res.time = calc_time_microseconds([&](){
		auto f = calc_right_part(u, c);
		boundary_setter_t set_boundary_conditions = bind(write_first_boundary_conditions, _1, _2, _3, _4, u);

		grid_t grid;
		grid.calc(x_grid, y_grid);

		vector_t q0 = calc_true_approx(bind(u, _1, _2, time_grid(0)), grid.bes);
		vector_t q1 = calc_true_approx(bind(u, _1, _2, time_grid(1)), grid.bes);
		vector_t q = calc_true_approx(bind(u, _1, _2, time_grid.back()), grid.bes);

		auto steps = solve_differential_equation(f, set_boundary_conditions, q0, q1, c, grid, time_grid);

		res.integral_residual = calc_integral_norm(bind(u, _1, _2, time_grid.back()), grid.es, steps.back());
		res.norm_residual = (q-steps.back()).norm() / q.size();
	});
	return res;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
template<class Ret, class Key>
class async_performer_t
{
public:
	void add(const function<Ret(void)>& f, const Key& key) {
		mf[key] = async(std::launch::deferred, f);
	}

	void finish(void) {
		int counter = 0;
		for (auto i = mf.begin(); i != mf.end(); ++i) {
			//if (counter % (mf.size()/10000 + 1) == 0)
				write_percent(double(counter)/mf.size());
			auto value = i->second.get();
			m[i->first] = value;
			counter++;
		}
		cout << "\r       \r";
	}

	auto begin(void) { return m.begin(); }
	auto end(void) { return m.end(); }

	Ret& operator[](const Key& key) { return m[key]; }
	const Ret& operator[](const Key& key) const { return m[key]; }
private:
	map<Key, future<Ret>> mf;
	map<Key, Ret> m;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
template<class ForwardIt, class GetValue>
double max_element_ignore_nan(ForwardIt first, ForwardIt last, GetValue get) {
	return get(*max_element(first, last, [get] (auto& a, auto& b) -> bool { 
		if (isnan(get(a)))
			return true;
		else
			return get(a) < get(b); 
	}));
}

//-----------------------------------------------------------------------------
void investigate_t_changing(
	int n,
	const string& filename,
	const function<fem_result_t(double)>& ft
) {
	auto uniform_value = ft(0);

	async_performer_t<fem_result_t, int> performer;

	grid_generator_t grid(-1, 1, n);
	for (int i = 0; i < grid.size(); i++) {
		performer.add([i, ft, grid] () -> fem_result_t {
			return ft(grid(i));
		}, i);
	}

	performer.finish();

	int counter = 0;
	auto integral_residual_max = max_element_ignore_nan(performer.begin(), performer.end(), [] (auto& a) -> double { return a.second.integral_residual; });
	auto norm_residual_max = max_element_ignore_nan(performer.begin(), performer.end(), [] (auto& a) -> double { return a.second.norm_residual; });
	
	ofstream fout(filename + ".txt");
	fout << "t\tintegral\tnorm\tuniform_integral\tuniform_norm\ttime" << endl;
	for (int i = 0; i < grid.size(); i++) {
		auto v = performer[i];
		fout 
			<< grid(i) << "\t"
			<< (isnan(v.integral_residual) ? integral_residual_max : v.integral_residual) << "\t"
			<< (isnan(v.norm_residual) ? norm_residual_max : v.norm_residual) << "\t"
			<< uniform_value.integral_residual << "\t"
			<< uniform_value.norm_residual << "\t"
			<< v.time << endl; 
	}
	fout.close();
}

//-----------------------------------------------------------------------------
void investigate_t2_changing(
	int n,
	const string& filename,
	const function<fem_result_t(double, double)>& ft
) {
	async_performer_t<fem_result_t, pair<int, int>> performer;

	grid_generator_t grid(-1, 1, n);
	for (int i = 0; i < grid.size(); i++) {
		for (int j = 0; j < grid.size(); j++) {
			performer.add([i, j, ft, grid] () -> fem_result_t {
				return ft(grid(i), grid(j));
			}, {i, j});
		}
	}

	performer.finish();

	auto integral_residual_max = max_element_ignore_nan(performer.begin(), performer.end(), [] (auto& a) -> double { return a.second.integral_residual; });
	auto norm_residual_max = max_element_ignore_nan(performer.begin(), performer.end(), [] (auto& a) -> double { return a.second.norm_residual; });

	ofstream fout(filename + ".integral.txt");
	ofstream fout2(filename + ".norm.txt");
	ofstream fout3(filename + ".time.txt");
	int last_line = 0;
	for (int i = 0; i < grid.size(); i++) {
		for (int j = 0; j < grid.size(); j++) {
			auto v = performer[{i, j}];
			fout << (isnan(v.integral_residual) ? integral_residual_max : v.integral_residual) << "\t";
			fout2 << (isnan(v.norm_residual) ? norm_residual_max : v.norm_residual) << "\t";
			fout3 << v.time << "\t";
		}
		fout << endl;
		fout2 << endl;
		fout3 << endl;
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

//-----------------------------------------------------------------------------
void investigate_functions(
	const string& filename,
	const function<fem_result_t(const function_3d_t&)>& f
) {
	vector<pair<function_3d_t, string>> spaces, times;

	spaces.push_back({[] (double x, double y, double t) -> double { return 1; }, "$1$"});
	spaces.push_back({[] (double x, double y, double t) -> double { return x+y; }, "$x+y$"});
	spaces.push_back({[] (double x, double y, double t) -> double { return x*x+y*y; }, "$x^2+y$"});
	spaces.push_back({[] (double x, double y, double t) -> double { return x*x*y+y*y*y; }, "$x^2y+y^3$"});
	spaces.push_back({[] (double x, double y, double t) -> double { return x*y*y; }, "$xy^2$"});
	spaces.push_back({[] (double x, double y, double t) -> double { return x*x*x*x+y*y*y*y; }, "$x^4+y^4$"});
	spaces.push_back({[] (double x, double y, double t) -> double { return exp(x*y); }, "$e^{xy}$"});

	times.push_back({[] (double x, double y, double t) -> double { return 0; }, "$0$"});
	times.push_back({[] (double x, double y, double t) -> double { return t; }, "$t$"});
	times.push_back({[] (double x, double y, double t) -> double { return t*t; }, "$t^2$"});
	times.push_back({[] (double x, double y, double t) -> double { return t*t*t; }, "$t^3$"});
	times.push_back({[] (double x, double y, double t) -> double { return t*t*t; }, "$t^4$"});
	times.push_back({[] (double x, double y, double t) -> double { return exp(t); }, "$e^t$"});

	async_performer_t<fem_result_t, pair<string, string>> performer;

	for (auto& i : spaces) {
		for (auto& j : times) {
			performer.add([i, j, &f]() -> fem_result_t {
				return f(function_3d_t([&] (double x, double y, double t) -> double { return i.first(x, y, t) + j.first(x, y, t); }));
			}, {i.second, j.second});
		}
	}

	performer.finish();
	
	ofstream fout(filename);
	fout << "a\t";
	for (auto& i : times)
		fout << i.second << "\t";
	for (auto& i : spaces) {
		fout << endl << i.second << "\t";
		for (auto& j : times) {
			auto v = performer[{i.second, j.second}];
			fout << "\\scalebox{.75}{\\tcell{$" << write_for_latex_double(v.integral_residual, 2) << "$\\\\$" << write_for_latex_double(v.norm_residual, 2) << "$\\\\$" << int(v.time/1000) << "$}}\t";
		}
	}
	fout.close();
}

//-----------------------------------------------------------------------------
void investigate_grid_changing(
	const string& filename,
	const function<pair<fem_result_t, int>(int)>& fi,
	int n
) {
	async_performer_t<pair<fem_result_t, int>, int> performer;

	for (int i = 0; i < n; i+=3) {
		performer.add([i, fi] () -> pair<fem_result_t, int> {
			return fi(i);
		}, i);
	}

	performer.finish();

	ofstream fout(filename + ".txt");
	fout << "i\tintegral\tnorm\ttime" << endl;
	for (int i = 0; i < n; i+=3) {
		auto v = performer[i];
		fout 
			<< v.second << "\t"
			<< v.first.integral_residual << "\t"
			<< v.first.norm_residual << "\t"
			<< v.first.time << endl; 
	}
	fout.close();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	cout << calc_time_microseconds([](){
		investigate_grid_changing(
			"space_sgrid",
			[] (int sz) -> pair<fem_result_t, int> {
				return {
					calc_fem_residual(
						[] (double x, double y, double t) -> double { return exp(x*y) + exp(t*t); }, 
						grid_generator_t(0, 1, 5+sz), 
						grid_generator_t(0, 1, 5+sz), 
						grid_generator_t(0, 1, 300)
					), 
					(7+sz)*(7+sz)
				};
			},
			70
		);

		investigate_grid_changing(
			"time_sgrid",
			[] (int sz) -> pair<fem_result_t, int> {
				return {
					calc_fem_residual(
						[] (double x, double y, double t) -> double { return exp(x*y) + exp(t*t); }, 
						grid_generator_t(0, 1, 20), 
						grid_generator_t(0, 1, 20), 
						grid_generator_t(0, 1, 1+sz)
					), 
					3+sz
				};
			},
			500
		);

		investigate_functions(
			"functions_table_10_10_10.txt",
			[] (const function_3d_t& u) -> fem_result_t {
				return calc_fem_residual(u, grid_generator_t(0, 1, 10), grid_generator_t(0, 1, 10), grid_generator_t(0, 1, 10));
			}
		);

		investigate_functions(
			"functions_table_50_50_50.txt",
			[](const function_3d_t& u) ->  fem_result_t {
				return calc_fem_residual(u, grid_generator_t(0, 1, 50), grid_generator_t(0, 1, 50), grid_generator_t(0, 1, 50));
			}
		);

		vector<pair<function_3d_t, int>> u_space_mas;
		u_space_mas.push_back({[] (double x, double y, double t) -> double { return x*x + y*y + t*t; }, 0});
		u_space_mas.push_back({[] (double x, double y, double t) -> double { return x*x*x*x + y*y*y*x + t*t*t*t; }, 1});
		u_space_mas.push_back({[] (double x, double y, double t) -> double { return exp(x*y) + exp(t*t); }, 2});
		u_space_mas.push_back({[] (double x, double y, double t) -> double { return exp((1-x)*(1-y)) + exp((1-t)*(1-t)); }, 3});
		u_space_mas.push_back({[] (double x, double y, double t) -> double { return x*x*x + y*y*y*y*x*x*t + t*t*exp(t); }, 4});
		u_space_mas.push_back({ [](double x, double y, double t) -> double { return x * x * x + y * y * y + t * t; }, 0 });

		for (auto& i : u_space_mas) {
			auto& u = i.first;
			investigate_t_changing(
				750,
				"time_tgrid_" + to_string(i.second),
				[u] (double tt) -> fem_result_t {
					return calc_fem_residual(u, grid_generator_t(0, 1, 10), grid_generator_t(0, 1, 10), grid_generator_t(0, 1, 10, tt));
				}
			);

			investigate_t2_changing(
				75,
				"space_tgrid_" + to_string(i.second),
				[u] (double tx, double ty) -> fem_result_t {
					return calc_fem_residual(u, grid_generator_t(0, 1, 10, tx), grid_generator_t(0, 1, 10, ty), grid_generator_t(0, 1, 10));
				}
			);
		}
	})/1000/1000 << "s" << endl;
	system("pause");
}