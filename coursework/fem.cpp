#include "fem.h"

//-----------------------------------------------------------------------------
bool basic_elem_t::is_boundary(void) const {
	return 
		up == nullptr ||
		down == nullptr ||
		left == nullptr ||
		right == nullptr;
}

//-----------------------------------------------------------------------------
double elem_t::get_hx(void) const { 
	return e[1]->x - e[0]->x; 
}

//-----------------------------------------------------------------------------
double elem_t::get_hy(void) const { 
	return e[3]->y - e[0]->y; 
}

//-----------------------------------------------------------------------------
double elem_t::value(double x, double y, const vector_t& q) const {
	// TODO
	return 0;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double non_linear_offset(double x, double t) {
	int signt = (t > 0) ? 1 : -1;
	t *= signt;
	t = 1.0 - t;
	t = (signt == -1) ? 1.0/t : t;
	return (1.0 - pow(t, x))/(1.0 - t);
}

//-----------------------------------------------------------------------------
grid_generator_t::grid_generator_t(double a, double b, int n, double t) : a(a), len(b-a), n1(n+1.0), t(t) {}

//-----------------------------------------------------------------------------
double grid_generator_t::operator()(int i) const {
	return a + len * non_linear_offset(i/n1, t);
}

//-----------------------------------------------------------------------------
int grid_generator_t::size() const {
	return n1+1;
}

//-----------------------------------------------------------------------------
void grid_t::calc(const grid_generator_t& gx, const grid_generator_t& gy) {
	bes.clear();
	bes.resize(gx.size() * gy.size());
	int counter = 0;
	for (int j = 0; j < gy.size(); ++j) {
		double y = gy(j);
		for (int i = 0; i < gx.size(); ++i) {
			double x = gx(i);
			basic_elem_t* down = (counter >= gx.size()) ? &bes[counter-gx.size()] : nullptr;
			basic_elem_t* left = (counter % gx.size() > 0) ? &bes[counter-1] : nullptr;
			bes[counter] = {counter, x, y, 
				nullptr, 
				down,
				left,
				nullptr
			};
			if (down != nullptr) down->up = &bes[counter];
			if (left != nullptr) left->right = &bes[counter];
			counter++;
		}
	}

	es.clear();
	counter = 0;
	for (auto& i : bes) {
		if (i.right != nullptr && i.up != nullptr && i.up->right == i.right->up && i.up->right != nullptr) {
			es.push_back({counter, 
				i.right->left, 
				i.right,
				i.up, 
				i.up->right
			});
			counter++;
		}
	}

	n = bes.size();
}

//-----------------------------------------------------------------------------
vector_t calc_true_approx(const function_2d_t& u, const vector<basic_elem_t>& bes) {
	vector_t result(bes.size());
	for (int i = 0; i < bes.size(); ++i)
		result(i) = u(bes[i].x, bes[i].y);
	return result;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
matrix_t calc_local_matrix_g(
	const elem_t& e, 
	const constants_t& cs
) {
	double hx = e.get_hx();
	double hy = e.get_hy();
	matrix_t result;
	matrix_t a(4, 4), b(4, 4);
	a <<
		2, -2, 1, -1,
		-2, 2, -1, 1,
		1, -1, 2, -2,
		-1, 1, -2, 2;
	b <<
		2, 1, -2, -1,
		1, 2, -1, -2,
		-2, -1, 2, 1,
		-1, -2, 1, 2;
	result = cs.lambda/6.0*(hy/hx*a + hx/hy*b);
	return result;
}

//-----------------------------------------------------------------------------
matrix_t calc_local_matrix_c(
	const elem_t& e
) {
	double hx = e.get_hx();
	double hy = e.get_hy();
	matrix_t result;
	matrix_t c(4, 4);
	c <<
		4, 2, 2, 1,
		2, 4, 1, 2,
		2, 1, 4, 2,
		1, 2, 2, 4;
	result = hx*hy/36.0*c;
	return result;
}

//-----------------------------------------------------------------------------
matrix_t calc_local_matrix_m(
	const elem_t& e,
	const constants_t& cs
) {
	return cs.gamma * calc_local_matrix_c(e);
}

//-----------------------------------------------------------------------------
vector_t calc_local_vector_b(
	const elem_t& e,
	const function_2d_t& f
) {
	vector_t fv(4);
	fv << 
		f(e.e[0]->x, e.e[0]->y),
		f(e.e[1]->x, e.e[1]->y),
		f(e.e[2]->x, e.e[2]->y),
		f(e.e[3]->x, e.e[3]->y);
	return calc_local_matrix_c(e) * fv;	
}

//-----------------------------------------------------------------------------
vector_t calc_global_vector(
	const vector<elem_t>& es,
	const function<vector_t(const elem_t&)> calc_local_vector,
	int n
) {
	vector_t result(n);
	result.fill(0);
	for (auto& e : es) {
		auto b = calc_local_vector(e);
		for (int i = 0; i < 4; ++i) {
			result(e.e[i]->i) += b(i);
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
matrix_sparse_t calc_global_matrix(
	const vector<elem_t>& es,
	const function<matrix_t(const elem_t&)> calc_local_matrix,
	int n
) {
	matrix_sparse_ra_t result(n);
	for (auto& e : es) {
		auto m = calc_local_matrix(e);
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				result(e.e[i]->i, e.e[j]->i) += m(i, j);
			}
		}
	}
	return result.to_sparse();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double calc_integral_gauss3(
	double a, double b, int n, // n - количество внутренных узлов
	const function_1d_t& f
) {
	const double x1 = -sqrt(3.0/5.0);
	const double x2 = 0;
	const double x3 = -x1;
	const double q1 = 5.0/9.0;
	const double q2 = 8.0/9.0;
	const double q3 = q1;
	double sum = 0;
	double xk = 0;
	double h = (b-a)/double(n+1);
	double h2 = h/2.0;
	
	for (int i = 0; i < n+1; ++i) {
		xk = a + h*i + h2;
		sum += q1 * f(xk + x1 * h2);
		sum += q2 * f(xk + x2 * h2);
		sum += q3 * f(xk + x3 * h2);
	}

	sum *= h;
	sum /= 2.0;
	return sum;
}

//-----------------------------------------------------------------------------
double calc_integral_gauss3(
	double ax, double bx, int nx,
	double ay, double by, int ny,
	const function_2d_t& f
) {
	return calc_integral_gauss3(ax, bx, nx, [ay, by, ny, f](double x)->double {
		return calc_integral_gauss3(ay, by, ny, bind(f, x, _1));
	});
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
function_1d_t calc_first_derivative(const function_1d_t& f) {
	return [f](double x) -> double {
		const double h = 0.001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}

//-----------------------------------------------------------------------------
function_1d_t calc_second_derivative(const function_1d_t& f) {
	return [f](double x) -> double {
		const double h = 0.001;
		return (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h))/(12*h*h);
	};
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
function_3d_t calc_right_part(
	const function_3d_t& u,
	const constants_t& cs
) {
	// f = -div(lambda * grad u) + gamma * u + sigma * du/dt + chi * d^2 u/dt^2
	return [=](double x, double y, double t) -> double {
		using namespace placeholders;
		auto ut = calc_first_derivative(bind(u, x, y, _1));

		auto uxx = calc_second_derivative(bind(u, _1, y, t));
		auto uyy = calc_second_derivative(bind(u, x, _1, t));
		auto utt = calc_second_derivative(bind(u, x, y, _1));

		return -cs.lambda * (uxx(x) + uyy(y)) + cs.gamma * u(x, y, t) + cs.sigma * ut(t) + cs.chi * utt(t);
	};
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void calc_crank_nicolson_method(
	const matrix_sparse_t& c,
	const matrix_sparse_t& m,
	const matrix_sparse_t& g,
	const vector_t& b0, // b current (b_j)
	const vector_t& bl, // b last (b_{j-1})
	const vector_t& bll, // b last last (b_{j-2})
	const vector_t& ql, // q last (q_{j-1})
	const vector_t& qll, // q last last (q_{j-2})
	const constants_t& cs,
	const grid_generator_t& time_grid,
	int time_i,
	matrix_sparse_t& a,
	vector_t& b
) {
	// Схема Кранка-Николсона
	//a = 0.5 * m + (cs.chi/dt/dt + cs.sigma/2.0/dt) * c + 0.5 * g;
	//b = 0.5 * (b0 + bll) - m * (cs.chi/dt/dt * (-2.0*ql + qll) - cs.sigma/2.0/dt * qll) - 0.5*g * qll - 0.5*m * qll; 

	// TODO

	// Моя простая схема
	//a = g  + m + (cs.chi/dt/dt + cs.sigma/2/dt) * c;
	//b = b0 + c * (2*cs.chi/dt/dt * ql + (cs.sigma/2/dt - cs.chi/dt/dt) * qll);
}

//-----------------------------------------------------------------------------
void write_first_boundary_conditions(
	matrix_sparse_t& a,
	vector_t& b,
	const vector<basic_elem_t>& bes,
	double t,
	const function_3d_t& u
) {
	for (int i = 0; i < bes.size(); ++i) {
		if (bes[i].is_boundary()) {
			a.clear_line(bes[i].i);
			a.d[bes[i].i] = 1;
			b(bes[i].i) = u(bes[i].x, bes[i].y, t);
		}
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
vector<vector_t> solve_differential_equation(
	const function_3d_t& f,
	const boundary_setter_t& set_boundary_conditions,
	const vector_t& q0,
	const vector_t& q1,
	const constants_t& cs,
	const grid_t& grid,
	const grid_generator_t& time_grid
) {
	auto c = calc_global_matrix(grid.es, calc_local_matrix_c, grid.n);
	auto m = calc_global_matrix(grid.es, bind(calc_local_matrix_m, _1, cs), grid.n);
	auto g = calc_global_matrix(grid.es, bind(calc_local_matrix_g, _1, cs), grid.n);

	auto calc_global_vector_b = [&] (int i) {
		return calc_global_vector(
			grid.es, 
			bind(
				calc_local_vector_b, 
				_1, 
				function_2d_t(bind(f, _1, _2, time_grid(i)))
			), 
			grid.n
		);
	};

	vector_t bll = calc_global_vector_b(0);
	vector_t bl = calc_global_vector_b(1);
	vector_t b0;

	vector_t qll = q0;
	vector_t ql = q1;
	vector_t q;

	vector<vector_t> result;
	result.push_back(q0);
	result.push_back(q1);

	matrix_sparse_t a(grid.n);
	vector_t b(grid.n);
	for (int i = 2; i < time_grid.size(); ++i) {
		b0 = calc_global_vector_b(i);
		calc_crank_nicolson_method(c, m, g, b0, bl, bll, ql, qll, cs, time_grid, i, a, b);
		set_boundary_conditions(a, b, grid.bes, time_grid(i));
		
		q = solve_by_los_lu(a, b, 1000, 1e-16, false);

		result.push_back(q);

		bll = bl;
		bl = b0;

		qll = ql;
		ql = q;
	}

	return result;
}
