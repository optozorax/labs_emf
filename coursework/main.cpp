#include <iostream>
#include <cmath>
#include "lib.h"
#include "fem.h"

using namespace std;
using namespace placeholders;

int main() {
	constants_t c = {1, 1, 1, 1};
	auto u = [] (double x, double y, double t) -> double { return exp(x*x*(y-t)) + t*t; };
	auto f = calc_right_part(u, c);
	boundary_setter_t set_boundary_conditions = bind(write_first_boundary_conditions, _1, _2, _3, _4, u);

	grid_t grid;
	grid.calc(grid_generator_t(0, 1, 1), grid_generator_t(0, 1, 1));

	grid_generator_t time_grid(0, 1, 1);

	vector_t q0 = calc_true_approx(bind(u, _1, _2, time_grid(0)), grid.bes);
	vector_t q1 = calc_true_approx(bind(u, _1, _2, time_grid(1)), grid.bes);
	vector_t q = calc_true_approx(bind(u, _1, _2, time_grid(time_grid.size()-1)), grid.bes);

	auto res = solve_differential_equation(f, set_boundary_conditions, q0, q1, c, grid, time_grid);

	double norm = (q-res.back()).norm();

	cout << "norm: " << norm << endl;
	cout << res.back() << endl;

	system("pause");
}