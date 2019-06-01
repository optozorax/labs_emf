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
double basic_function(double left, double middle, double right, double x) {
	// Базовая линейная функция
	if (x > left && x < right) {
		if (x < middle) {
			if (middle != left)
				return (x-left)/(middle-left);
			else
				return 0;
		} else {
			if (right != middle)
				return 1-(x-middle)/(right-middle);
			else
				return 0;
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
double basic_function_grad(double left, double middle, double right, double x) {
	// Базовая линейная функция
	if (x > left && x < right) {
		if (x < middle) {
			if (middle != left)
				return 1.0/(middle-left);
			else
				return 0;
		} else {
			if (right != middle)
				return -1.0/(right-middle);
			else
				return 0;
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

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
matrix_t calc_global_vector(
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
matrix_t calc_global_matrix(
	const vector<elem_t>& es,
	const function<matrix_t(const elem_t&)> calc_local_matrix,
	int n
) {
	matrix_t result(n, n);
	result.fill(0);
	for (auto& e : es) {
		auto m = calc_local_matrix(e);
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				result(e.e[i]->i, e.e[j]->i) += m(i, j);
			}
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
void calc_crank_nicolson_method(
	const matrix_t& c,
	const matrix_t& m,
	const matrix_t& g,
	const vector_t& b0, // b current (b_j)
	const vector_t& bl, // b last (b_{j-1})
	const vector_t& bll, // b last last (b_{j-2})
	const vector_t& ql, // q last (q_{j-1})
	const vector_t& qll, // q last last (q_{j-2})
	const constants_t& cs,
	double dt,
	matrix_t& a,
	vector_t& b
) {
	// Схема Кранка-Николсона
	a = 0.5 * m + (cs.chi/dt/dt + cs.sigma/2.0/dt) * c + 0.5 * g;
	b = 0.5 * (b0 + bll) - m * (cs.chi/dt/dt * (-2.0*ql + qll) - cs.sigma/2.0/dt * qll) - 0.5*g * qll - 0.5*m * qll;

	// Моя простая схема
	//a = g  + m + (cs.chi/dt/dt + cs.sigma/2/dt) * c;
	//b = b0 + c * (2*cs.chi/dt/dt * ql + (cs.sigma/2/dt - cs.chi/dt/dt) * qll);
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
		const double h = 0.001;
		return (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h))/(12*h*h);
	};
}

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
	const constants_t& cs,
	const grid_t& grid,
	int nt, double dt, double at
) {
	auto c = calc_global_matrix(grid.es, calc_local_matrix_c, grid.bn);
	auto m = calc_global_matrix(grid.es, bind(calc_local_matrix_m, _1, cs), grid.bn);
	auto g = calc_global_matrix(grid.es, bind(calc_local_matrix_g, _1, cs), grid.bn);

	vector_t bll = calc_global_vector(grid.es, bind(calc_local_vector_b, _1, function_2d_t(bind(f, _1, _2, at))), grid.bn);
	vector_t bl = calc_global_vector(grid.es, bind(calc_local_vector_b, _1, function_2d_t(bind(f, _1, _2, at+dt))), grid.bn);
	vector_t b0;

	vector_t qll = q0;
	vector_t ql = q1;
	vector_t q;

	vector<vector_t> result;
	result.push_back(q0);
	result.push_back(q1);

	matrix_t a(grid.bn, grid.bn);
	vector_t b(grid.bn);
	for (int i = 2; i < nt; ++i) {
		double t = at + i*dt;
		b0 = calc_global_vector(grid.es, bind(calc_local_vector_b, _1, function_2d_t(bind(f, _1, _2, t))), grid.bn);
		calc_crank_nicolson_method(c, m, g, b0, bl, bll, ql, qll, cs, dt, a, b);
		set_boundary_conditions(a, b, grid.bes, t);
		q = Eigen::JacobiSVD<matrix_t>(a, Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

		//cout << a << endl;
		//cout << b << endl;

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
	auto u = [] (double x, double y, double t) -> double { return exp(x*x*(y-t)) + t*t; };
	auto f = calc_right_part(u, c);
	boundary_setter_t set_boundary_conditions = bind(write_first_boundary_conditions, _1, _2, _3, _4, u);

	grid_t grid;
	grid.calc(0, 0, 1, 1, 10, 10);

	double at = 0, bt = 1, nt = 10;
	double dt = (bt-at)/(nt-1);

	vector_t q0 = calc_true_approx(bind(u, _1, _2, at), grid.bes);
	vector_t q1 = calc_true_approx(bind(u, _1, _2, at+dt), grid.bes);
	vector_t q = calc_true_approx(bind(u, _1, _2, bt), grid.bes);

	auto res = solve_differential_equation(f, set_boundary_conditions, q0, q1, c, grid, nt, dt, at);

	double norm = (q-res.back()).norm();

	cout << "norm: " << norm << endl;
	cout << res.back() << endl;

	system("pause");
}