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
	return e[2]->x - e[0]->x; 
}

//-----------------------------------------------------------------------------
double elem_t::get_hy(void) const { 
	return e[6]->y - e[0]->y; 
}

//-----------------------------------------------------------------------------
double elem_t::value(double x, double y, const vector_t& q) const {
	double xp = e[0]->x;
	double yp = e[0]->y;
	double hx = get_hx();
	double hy = get_hy();
	auto x1 = [&](double x) -> double { return 1.0 + 3.0*xp/hx + 2.0*(x*x)/(hx*hx) - x*(3.0*hx + 4.0*xp)/(hx*hx) + 2.0*(xp*xp)/(hx*hx); };
	auto x2 = [&](double x) -> double { return -4.0*(x*x)/(hx*hx) + 4.0*x*(hx + 2.0*xp)/(hx*hx) - 4.0*xp*(hx + xp)/(hx*hx); };
	auto x3 = [&](double x) -> double { return 2.0*(x*x)/(hx*hx) - x*(hx + 4.0*xp)/(hx*hx) + xp*(hx + 2.0*xp)/(hx*hx); };
	auto y1 = [&](double y) -> double { return 1.0 + 3.0*yp/hy + 2.0*(y*y)/(hy*hy) - y*(3.0*hy + 4.0*yp)/(hy*hy) + 2.0*(yp*yp)/(hy*hy); };
	auto y2 = [&](double y) -> double { return -4.0*(y*y)/(hy*hy) + 4.0*y*(hy + 2.0*yp)/(hy*hy) - 4.0*yp*(hy + yp)/(hy*hy); };
	auto y3 = [&](double y) -> double { return 2.0*(y*y)/(hy*hy) - y*(hy + 4.0*yp)/(hy*hy) + yp*(hy + 2.0*yp)/(hy*hy); };

	auto psi1 = [&]() -> double { return x1(x) * y1(y); };
	auto psi2 = [&]() -> double { return x2(x) * y1(y); };
	auto psi3 = [&]() -> double { return x3(x) * y1(y); };
	auto psi4 = [&]() -> double { return x1(x) * y2(y); };
	auto psi5 = [&]() -> double { return x2(x) * y2(y); };
	auto psi6 = [&]() -> double { return x3(x) * y2(y); };
	auto psi7 = [&]() -> double { return x1(x) * y3(y); };
	auto psi8 = [&]() -> double { return x2(x) * y3(y); };
	auto psi9 = [&]() -> double { return x3(x) * y3(y); };

	double v1 = psi1() * q[e[0]->i];
	double v2 = psi2() * q[e[1]->i];
	double v3 = psi3() * q[e[2]->i];
	double v4 = psi4() * q[e[3]->i];
	double v5 = psi5() * q[e[4]->i];
	double v6 = psi6() * q[e[5]->i];
	double v7 = psi7() * q[e[6]->i];
	double v8 = psi8() * q[e[7]->i];
	double v9 = psi9() * q[e[8]->i];

	return v1 + v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9;
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
	if (t == 1.0) return x;
	return (1.0 - pow(t, x))/(1.0 - t);
}

//-----------------------------------------------------------------------------
grid_generator_t::grid_generator_t(double a, double b, int n, double t) : a(a), len(b-a), n1(2 * n+1.0), t(t) {}

//-----------------------------------------------------------------------------
double grid_generator_t::operator()(int i) const {
	return a + len * non_linear_offset(i/n1, t);
}

//-----------------------------------------------------------------------------
int grid_generator_t::size(void) const {
	return n1+1;
}

//-----------------------------------------------------------------------------
double grid_generator_t::back(void) const {
	return operator()(size()-1);
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
				nullptr,
			false
			};
			if (down != nullptr) down->up = &bes[counter];
			if (left != nullptr) left->right = &bes[counter];
			counter++;
		}
	}

	n = bes.size();

	es.clear();
	counter = 0;
	for (auto& i : bes) {
		bool is_has_next_elems = 
			i.right != nullptr && 
			i.right->right != nullptr && 
			i.up != nullptr && 
			i.up->up != nullptr && 
			i.up->up->right->right != nullptr;
		bool is_not_part_of_elem = i.is_part_of_elem != true;
		if (is_has_next_elems && is_not_part_of_elem) {
			es.push_back({counter, 
				i.right->left, 
				i.right,
				i.right->right, 
				i.up,
				i.up->right,
				i.up->right->right,
				i.up->up,
				i.up->up->right,
				i.up->up->right->right,
			});
			i.right->left->is_part_of_elem = true;
			i.right->is_part_of_elem = true;
			//i.right->right->is_part_of_elem = true;
			i.up->is_part_of_elem = true;
			i.up->right->is_part_of_elem = true;
			i.up->right->right->is_part_of_elem = true;
			//i.up->up->is_part_of_elem = true;
			i.up->up->right->is_part_of_elem = true;
			//i.up->up->right->right->is_part_of_elem = true;
			counter++;
		}
	}
}

//-----------------------------------------------------------------------------
vector_t calc_true_approx(const function_2d_t& u, const vector<basic_elem_t>& bes) {
	vector_t result(bes.size());
	for (int i = 0; i < bes.size(); ++i)
		result(i) = u(bes[i].x, bes[i].y);
	return result;
}

//-----------------------------------------------------------------------------
double calc_integral_norm(const function_2d_t& u, const vector<elem_t>& es, const vector_t& q) {
	double sum = 0;
	for (auto& i : es) {
		sum += calc_integral_gauss3(
			i.e[0]->x, i.e[2]->x, 10, 
			i.e[0]->y, i.e[6]->y, 10,
			[&](double x, double y) -> double {
				return fabs(u(x, y) - i.value(x, y, q));
			}
		);
	}
	return sum;
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
	matrix_t a(9, 9), b(9, 9);
	a <<
		28, -32, 4, 14, -16, 2, -7, 8, -1,
		-32, 64, -32, -16, 32, -16, 8, -16, 8,
		4, -32, 28, 2, -16, 14, -1, 8, -7,
		14, -16, 2, 112, -128, 16, 14, -16, 2,
		-16, 32, -16, -128, 256, -128, -16, 32, -16,
		2, -16, 14, 16, -128, 112, 2, -16, 14,
		-7, 8, -1, 14, -16, 2, 28, -32, 4,
		8, -16, 8, -16, 32, -16, -32, 64, -32,
		-1, 8, -7, 2, -16, 14, 4, -32, 28;
	b <<
		28, 14, -7, -32, -16, 8, 4, 2, -1,
		14, 112, 14, -16, -128, -16, 2, 16, 2,
		-7, 14, 28, 8, -16, -32, -1, 2, 4,
		-32, -16, 8, 64, 32, -16, -32, -16, 8,
		-16, -128, -16, 32, 256, 32, -16, -128, -16,
		8, -16, -32, -16, 32, 64, 8, -16, -32,
		4, 2, -1, -32, -16, 8, 28, 14, -7,
		2, 16, 2, -16, -128, -16, 14, 112, 14,
		-1, 2, 4, 8, -16, -32, -7, 14, 28;
	result = cs.lambda/90.0*(hy/hx*a + hx/hy*b);
	return result;
}

//-----------------------------------------------------------------------------
matrix_t calc_local_matrix_c(
	const elem_t& e
) {
	double hx = e.get_hx();
	double hy = e.get_hy();
	matrix_t result;
	matrix_t c(9, 9);
	c <<
		16, 8, -4, 8, 4, -2, -4, -2, 1,
		8, 64, 8, 4, 32, 4, -2, -16, -2,
		-4, 8, 16, -2, 4, 8, 1, -2, -4,
		8, 4, -2, 64, 32, -16, 8, 4, -2,
		4, 32, 4, 32, 256, 32, 4, 32, 4,
		-2, 4, 8, -16, 32, 64, -2, 4, 8,
		-4, -2, 1, 8, 4, -2, 16, 8, -4,
		-2, -16, -2, 4, 32, 4, 8, 64, 8,
		1, -2, -4, -2, 4, 8, -4, 8, 16;
	result = hx*hy/900.0*c;
	return result;
}

//-----------------------------------------------------------------------------
vector_t calc_local_vector_b(
	const elem_t& e,
	const function_2d_t& f
) {
	vector_t fv(9);
	fv << 
		f(e.e[0]->x, e.e[0]->y),
		f(e.e[1]->x, e.e[1]->y),
		f(e.e[2]->x, e.e[2]->y),

		f(e.e[3]->x, e.e[3]->y),
		f(e.e[4]->x, e.e[4]->y),
		f(e.e[5]->x, e.e[5]->y),

		f(e.e[6]->x, e.e[6]->y),
		f(e.e[7]->x, e.e[7]->y),
		f(e.e[8]->x, e.e[8]->y);
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
		for (int i = 0; i < 9; ++i) {
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
		for (int i = 0; i < 9; ++i) {
			for (int j = 0; j < 9; ++j) {
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
	const matrix_sparse_t& g,
	const vector_t& b0,
	const vector_t& bl,
	const vector_t& bll,
	const vector_t& ql,
	const vector_t& qll,
	const constants_t& cs,
	const grid_generator_t& time_grid,
	int time_i,
	matrix_sparse_t& a,
	vector_t& b
) {
	// Схема Кранка-Николсона

	// Константы для вычислений с неравномерной сеткой по времени
	double t0 = time_grid(time_i);
	double t1 = time_grid(time_i-1);
	double t2 = time_grid(time_i-2);

	double d1 = t0-t2;
	double d2 = (t0*(t0-t2-t1)+t2*t1)/2.0;
	double m1 = (t0-t2)/(t1-t2);
	double m2 = (t0-t1)/(t1-t2);

	// Вычисляемт матрицу a
	a = c;
	double c1 = cs.gamma/2.0 + cs.sigma/d1 + cs.chi/d2;
	for (int i = 0; i < a.d.size(); i++)
		a.d[i] = g.d[i]/2.0 + c.d[i]*c1;
	for (int i = 0; i < a.l.size(); i++) {
		a.l[i] = g.l[i]/2.0 + c.l[i]*c1;
		a.u[i] = g.u[i]/2.0 + c.u[i]*c1;
	}

	// Рассчитываем вектор b
	b = (b0 + bll)/2.0;
	vector_t temp = qll;
	g.mul(temp);
	b = b - temp/2.0;

	temp = ql*(m1*cs.chi/d2) + qll*(-cs.gamma/2.0 + cs.sigma/d1 - m2*cs.chi/d2);
	c.mul(temp);

	b = b + temp;
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
		calc_crank_nicolson_method(c, g, b0, bl, bll, ql, qll, cs, time_grid, i, a, b);
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
