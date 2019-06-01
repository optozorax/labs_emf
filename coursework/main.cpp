#include <iostream>
#include <functional>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace placeholders;

typedef Eigen::MatrixXd matrix_t;
typedef Eigen::VectorXd vector_t;

typedef function<double(double)> function_1d_t;
typedef function<double(double, double)> function_2d_t;
typedef function<double(double, double, double)> function_3d_t;

/** Вес, домноженный на базовую функцию. Из сумм этих элементов образуется конечный элемент. */
struct basic_elem_t
{
	int i;

	double x, y;

	basic_elem_t *up, *down, *left, *right;
};

bool is_boundary(const basic_elem_t& be) {
	return 
		be.up == nullptr ||
		be.down == nullptr ||
		be.left == nullptr ||
		be.right == nullptr;
}

/** Прямоугольный конечный элемент. */
struct elem_t
{
	int i;
	basic_elem_t* e[4];

	double get_hx(void) const { return e[1]->x - e[0]->x; }
	double get_hy(void) const { return e[3]->y - e[0]->y; }
};

/** Все константы решаемого уравнения. */
struct constants_t
{
	double lambda;
	double gamma;
	double sigma;
	double chi;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
matrix_t calc_local_matrix_g(
	const elem_t& e, 
	const constants_t& c
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
	result = c.lambda/6.0*(hy/hx*a + hx/hy*b);
	return result;
}

//-----------------------------------------------------------------------------
matrix_t calc_local_matrix_m(
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
vector_t calc_local_vector_b(
	const elem_t& e,
	const function_2d_t& f
) {
	double hx = e.get_hx();
	double hy = e.get_hy();
	matrix_t c(4, 4);
	c <<
		4, 2, 2, 1,
		2, 4, 1, 2,
		2, 1, 4, 2,
		1, 2, 2, 4;
	vector_t fv(4);
	fv << 
		f(e.e[0]->x, e.e[0]->y),
		f(e.e[1]->x, e.e[1]->y),
		f(e.e[2]->x, e.e[2]->y),
		f(e.e[3]->x, e.e[3]->y);
	vector_t result = hx*hy/36.0*c*fv;
	return result;	
}

//-----------------------------------------------------------------------------
matrix_t calc_global_matrix_g(
	const vector<elem_t>& es, 
	const constants_t& c,
	int n
) {
	matrix_t result(n, n);
	for (auto& e : es) {
		auto g = calc_local_matrix_g(e, c);
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				result(e.e[i]->i, e.e[j]->i) = g(i, j);
			}
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
matrix_t calc_global_matrix_m(
	const vector<elem_t>& es,
	int n
) {
	matrix_t result(n, n);
	for (auto& e : es) {
		auto m = calc_local_matrix_m(e);
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				result(e.e[i]->i, e.e[j]->i) = m(i, j);
			}
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
vector_t calc_global_vector_b(
	const vector<elem_t>& es,
	const function_2d_t& f,
	int n
) {
	vector_t result(n);
	for (auto& e : es) {
		auto b = calc_local_vector_b(e, f);
		for (int i = 0; i < 4; ++i) {
			result(e.e[i]->i) = b(i);
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
void calc_crank_nicolson_method(
	const matrix_t& m,
	const matrix_t& g,
	const vector_t& b0, // b current (b_j)
	const vector_t& bl, // b last (b_{j-1})
	const vector_t& bll, // b last last (b_{j-2})
	const vector_t& ql, // q last (q_{j-1})
	const vector_t& qll, // q last last (q_{j-2})
	const constants_t& c,
	double dt,
	matrix_t& a,
	vector_t& b
) {
	a = (c.chi/dt/dt + c.sigma/2.0/dt) * m + 0.5 * g;
	b = 0.5 * (b0 + bll) - m * (c.chi/dt/dt * (-2.0*ql + qll) - c.sigma/2.0/dt * qll) - 0.5*g * qll;

	//a = m * ((c.chi / dt / dt) + c.sigma / 2 / dt);
	//b = bl - m * (-2.0*c.chi/dt/dt * ql + c.chi/dt/dt * qll - c.sigma/2.0/dt * qll) - g*ql;
}

//-----------------------------------------------------------------------------
void write_first_boundary_conditions(
	matrix_t& a,
	vector_t& b,
	const vector<basic_elem_t>& bes,
	double t,
	const function_3d_t& u
) {
	auto clear_line = [&a] (int line) {
		for (int i = 0; i < a.cols(); ++i) {
			a(line, i) = 0;
		}
	};
	for (int i = 0; i < bes.size(); ++i) {
		if (is_boundary(bes[i])) {
			clear_line(bes[i].i);
			a(bes[i].i, bes[i].i) = 1;
			b(bes[i].i) = u(bes[i].x, bes[i].y, t);
		}
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

class grid_t
{
public:
	vector<elem_t> es;
	vector<basic_elem_t> bes;
	int bn;
	int n;

	void calc(double x0, double y0, double x1, double y1, int nx, int ny) {
		bes.clear();
		bes.resize(nx * ny);
		int counter = 0;
		for (int j = 0; j < ny; ++j) {
			double y = y0 + (y1-y0)/(ny-1)*j;
			for (int i = 0; i < nx; ++i) {
				double x = x0 + (x1-x0)/(nx-1)*i;
				basic_elem_t* down = (counter >= nx) ? &bes[counter-nx] : nullptr;
				basic_elem_t* left = (counter%nx > 0) ? &bes[counter-1] : nullptr;
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

		bn = bes.size();
		n = es.size();
	}
};

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
		const double h = 0.0001;
		return (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h))/(12*h*h);
	};
}

//-----------------------------------------------------------------------------
function_3d_t calc_right_part(
	const function_3d_t& u,
	const constants_t& c
) {
	// f = -div(lambda * grad u) + gamma * u + sigma * du/dt + chi * d^2 u/dt^2
	return [=](double x, double y, double t) -> double {
		using namespace placeholders;
		auto ut = calc_first_derivative(bind(u, x, y, _1));

		auto uxx = calc_second_derivative(bind(u, _1, y, t));
		auto uyy = calc_second_derivative(bind(u, x, _1, t));
		auto utt = calc_second_derivative(bind(u, x, y, _1));

		return -c.lambda * (uxx(x) + uyy(y)) + c.gamma * u(x, y, t) + c.sigma * ut(t) + c.chi * utt(t);
	};
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
typedef function<void(matrix_t&, vector_t&, const vector<basic_elem_t>&, double)> boundary_setter_t;

//-----------------------------------------------------------------------------
vector<vector_t> solve_differential_equation(
	const function_3d_t& f,
	const boundary_setter_t& set_boundary_conditions,
	const vector_t& q0,
	const vector_t& q1,
	const constants_t& c,
	const grid_t& grid,
	int nt, double dt, double at
) {
	auto m = calc_global_matrix_m(grid.es, grid.bn);
	auto g = calc_global_matrix_g(grid.es, c, grid.bn);

	vector_t bll = calc_global_vector_b(grid.es, bind(f, _1, _2, at), grid.bn);
	vector_t bl = calc_global_vector_b(grid.es, bind(f, _1, _2, at+dt), grid.bn);
	vector_t b0;

	vector_t qll = q0;
	vector_t ql = q1;
	vector_t q;

	vector<vector_t> result;
	result.push_back(q0);
	result.push_back(q1);

	matrix_t a;
	vector_t b;
	for (int i = 2; i < nt; ++i) {
		double t = at + i*dt;
		b0 = calc_global_vector_b(grid.es, bind(f, _1, _2, t), grid.bn);
		calc_crank_nicolson_method(m, g, b0, bl, bll, ql, qll, c, dt, a, b);
		set_boundary_conditions(a, b, grid.bes, t);
		q = Eigen::JacobiSVD<matrix_t>(a, Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

		result.push_back(q);

		bll = bl;
		bl = b0;

		qll = ql;
		ql = q;
	}

	return result;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	constants_t c = {1, 1, 1, 1};
	auto u = [] (double x, double y, double t) -> double { return x + y + t; };
	auto f = calc_right_part(u, c);
	boundary_setter_t set_boundary_conditions = bind(write_first_boundary_conditions, _1, _2, _3, _4, u);

	grid_t grid;
	grid.calc(0, 0, 1, 1, 3, 3);

	double at = 0, bt = 1, nt = 3;
	double dt = (bt-at)/(nt-1);

	vector_t q0 = calc_true_approx(bind(u, _1, _2, at), grid.bes);
	vector_t q1 = calc_true_approx(bind(u, _1, _2, at+dt), grid.bes);
	vector_t q = calc_true_approx(bind(u, _1, _2, bt), grid.bes);

	auto res = solve_differential_equation(f, set_boundary_conditions, q0, q1, c, grid, nt, dt, at);

	double norm = (q-res.back()).norm();

	cout << "norm: " << norm << endl;

	system("pause");
}