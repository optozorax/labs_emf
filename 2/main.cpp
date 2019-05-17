#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

typedef std::function<double(double)> Function1D;
typedef std::function<double(double, double)> Function2D;
typedef std::function<double(double, double, double)> Function3D;

//-----------------------------------------------------------------------------
double integral_gauss3(const std::vector<double>& X, Function1D f) {
	const double x1 = -sqrt(3.0/5.0);
	const double x2 = 0;
	const double x3 = -x1;
	const double q1 = 5.0/9.0;
	const double q2 = 8.0/9.0;
	const double q3 = q1;
	double sum = 0;
	double xk = 0;
	double h = X[1] - X[0];
	double h2 = h/2.0;
	
	for (int i = 0; i < X.size()-1; ++i) {
		xk = (X[i]+X[i+1])/2.0;
		sum += q1 * f(xk + x1 * h2);
		sum += q2 * f(xk + x2 * h2);
		sum += q3 * f(xk + x3 * h2);
	}

	sum *= h;
	sum /= 2.0;
	return sum;
}

//-----------------------------------------------------------------------------
void make_grid(std::vector<double>& X, double a, double b, long n) {
	X.clear();
	double h = (b-a)/n;
	for (double i = 0; i <= n; ++i)
		X.push_back(a + h * i);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double basicFunction(double left, double middle, double right, double x) {
	// Базовая линейная функция
	if (x > left && x < right) {
		if (x < middle) {
			return (x-left)/(middle-left);
		} else {
			return 1-(x-middle)/(right-middle);
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
double basicFunctionGrad(double left, double middle, double right, double x) {
	// Базовая линейная функция
	if (x > left && x < right) {
		if (x < middle) {
			return 1.0/(middle-left);
		} else {
			return -1.0/(right-middle);
		}
	}

	return 0;
}

class lin_approx_t
{
public:
	vector<double> q; /// Вектор весов
	vector<double> x; /// Массив положений каждого элемента

	int size(void) const {
		return x.size();
	}

	double left(int i) const {
		if (i == 0)
			return x[0]-0.00001;
		else
			return x[i-1];
	}

	double middle(int i) const {
		return x[i];
	}

	double right(int i) const {
		if (i == x.size() - 1)
			return x.back()+0.00001;
		else
			return x[i+1];
	}

	double value(double pos) const {
		int start = std::distance(x.begin(), std::lower_bound(x.begin(), x.end(), pos))-1;
		if (start == -1) {
			if (std::fabs(pos-x[0]) < 0.00001)
				return q[0];
			else
				return 0;
		}
		if (start == x.size()) {
			if (std::fabs(pos-x.back()) < 0.00001)
				return q.back();
			else
				return 0;
		}
		double sum = 0;
		for (int i = start; i < std::min<int>(start+2, x.size()); ++i)
		//for (int i = 0; i < x.size(); ++i)
			sum += q[i] * basic(pos, i);
		return sum;
	} /// Получить значение функции аппроксимации в некоторой точке

	double basic(double pos, int i) const {
		return basicFunction(left(i), middle(i), right(i), pos);
	}

	double basic_grad(double pos, int i) const {
		return basicFunctionGrad(left(i), middle(i), right(i), pos);	
	}
};

double calc_a1_integral(const Function1D& lambda, const lin_approx_t& u, int i, int j) {
	if (std::abs(i - j) > 1)
		return 0;
	auto f = [&] (double x) -> double {
		return lambda(u.value(x)) * u.basic_grad(x, i) * u.basic_grad(x, j);
	};
	std::vector<double> X;
	if (i != j)
		make_grid(X, std::min(u.middle(i), u.middle(j)), std::min(u.right(i), u.right(j)), 10);
	else
		make_grid(X, u.left(i), u.right(i), 10);
	return integral_gauss3(X, f);
}

double calc_a2_integral(const Function1D& sigma, const lin_approx_t& u, int i, int j) {
	if (std::abs(i - j) > 1)
		return 0;
	auto f = [&] (double x) -> double {
		return sigma(x) * u.basic(x, i) * u.basic(x, j);
	};
	std::vector<double> X;
	if (i != j)
		make_grid(X, std::min(u.middle(i), u.middle(j)), std::min(u.right(i), u.right(j)), 10);
	else
		make_grid(X, u.left(i), u.right(i), 10);
	return integral_gauss3(X, f);	
}

double calc_b1_integral(const Function1D& f, const lin_approx_t& u, int i) {
	auto fun = [&] (double x) -> double {
		return f(x) * u.basic(x, i);
	};
	std::vector<double> X;
	make_grid(X, u.left(i), u.right(i), 10);
	return integral_gauss3(X, fun);	
}

double calc_b2_integral(const Function1D& sigma, const lin_approx_t& u, const lin_approx_t& u_last_time, int i) {
	auto fun = [&] (double x) -> double {
		return sigma(x) * u_last_time.value(x) * u.basic(x, i);
	};
	std::vector<double> X;
	make_grid(X, u.left(i), u.right(i), 10);
	return integral_gauss3(X, fun);	
}

lin_approx_t approximate_function(Function1D lambda, const lin_approx_t& u) {
	lin_approx_t result = u;
	for (int i = 0; i < result.q.size(); ++i)
		result.q[i] = lambda(u.x[i]);
	return result;
}

Matrix calcA(const lin_approx_t& u_last, double dt, Function1D lambda, Function1D sigma) {
	Matrix result(u_last.size(), u_last.size());
	for (int i = 0; i < u_last.size()-1; ++i) {
		for (int j = 0; j < u_last.size(); ++j) {
			result(i, j) = calc_a1_integral(lambda, u_last, i, j) + calc_a2_integral(sigma, u_last, i, j) / dt;
		}
	}

	return result;
} /// Рассчитать матрицу A(q)
Vector calcB(const lin_approx_t& u_last, double dt, Function1D f, Function1D sigma, const lin_approx_t& u_last_time) {
	Vector result(u_last.size());
	for (int i = 0; i < u_last.size(); ++i) {
		result(i) = calc_b1_integral(f, u_last, i) + calc_b2_integral(sigma, u_last, u_last_time, i) / dt;
	}

	return result;
} /// Рассчитать вектор правой части b(q)

void write_first_boundary_conditions(Matrix& A, Vector& b, const lin_approx_t& u, Function2D u_true, double time) {
	A(0, 0) = 1;
	for (int i = 1; i < A.cols(); ++i) A(0, i) = 0;
	b(0) = u_true(u.x.front(), time);

	A(A.rows()-1, A.cols()-1) = 1;
	for (int i = 0; i < A.cols() - 1; ++i) A(A.rows()-1, i) = 0;
	b(b.rows()-1) = u_true(u.x.back(), time);
}

struct Result
{
	enum ExitType {
		EXIT_RESIDUAL,
		EXIT_STEP,
		EXIT_ITERATIONS,
		EXIT_ERROR
	} exitType;

	lin_approx_t answer;
	int iterations;
	double residual;
};

ostream& operator<<(ostream& out, const Result::ExitType& exit) {
	switch (exit) {
		case Result::EXIT_RESIDUAL: out << "residual"; break;
		case Result::EXIT_STEP: out << "step"; break;
		case Result::EXIT_ITERATIONS: out << "iterations"; break;
		case Result::EXIT_ERROR: out << "error"; break;
	}
	return out;
}

Result solveFixedPointIteration(
	Function2D f,
	Function2D u_true,
	Function1D lambda,
	Function1D sigma,
	const lin_approx_t& u_last_time, 
	double dt,
	double time,
	double eps, 
	int maxiter
) {
	lin_approx_t u = u_last_time;

	Vector last_x(u.q.size());
	for (int i = 0; i < u.q.size(); i++)
		last_x(i) = u.q[i];

	Result result;
	result.iterations = 0;
	while (true) {
		auto A = calcA(u, dt, lambda, sigma);
		auto b = calcB(u, dt, std::bind(f, std::placeholders::_1, time), sigma, u_last_time);
		write_first_boundary_conditions(A, b, u, u_true, time);

		//cout << "x:" << endl << last_x << endl;
		//cout << "A: " << endl << A << endl;
		//cout << "b: " << endl << b << endl;

		Eigen::JacobiSVD<Matrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Vector x = svd.solve(b);
		for (int i = 0; i < u.size(); ++i)
			u.q[i] = x(i);

		result.iterations++;

		if (result.iterations > maxiter) {
			result.exitType = Result::EXIT_ITERATIONS;
			break;
		}
		
		{
			auto A1 = calcA(u, dt, lambda, sigma);
			auto b1 = calcB(u, dt, std::bind(f, std::placeholders::_1, time), sigma, u_last_time);
			write_first_boundary_conditions(A1, b1, u, u_true, time);

			double au = (A1 * x - b1).norm();
			double bu = b1.norm();
			result.residual = au/bu;
			if (au/bu < eps) {
				result.exitType = Result::EXIT_RESIDUAL;
				break;
			}
		}

		if ((x-last_x).norm() < eps) {
			result.exitType = Result::EXIT_STEP;
			break;
		}

		last_x = x;
	}

	result.answer = u;

	return result;
} /// Решить уравнение методом простой итерации

Function1D calcFirstDerivative(const Function1D& f) {
	return [f](double x) -> double {
		const double h = 0.001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}

Function2D calcRightPart(
	const Function1D& lambda, 
	const Function2D& u, 
	const Function1D& sigma
) {
	// f = -div(lambda(u) * grad u) + sigma * du/dt
	return [=](double x, double t) -> double {
		using namespace std::placeholders;
		auto ut = calcFirstDerivative(std::bind(u, x, _1));
		auto ux = calcFirstDerivative(std::bind(u, _1, t));
		auto lambda_grad = [lambda, ux, u](double x, double t) -> double {
			return lambda(u(x, t)) * ux(x);
		};
		auto div = calcFirstDerivative(std::bind(lambda_grad, _1, t));
		return -div(x) + sigma(x) * ut(x);
	};
}

double norm(const std::vector<double>& a, const std::vector<double>& b) {
	assert(a.size() == b.size());
	double sum = 0;
	for (int i = 0; i < a.size(); i++)
		sum += (a[i] - b[i])*(a[i] - b[i]);
	return sqrt(sum);
}

ostream& operator<<(ostream& out, const std::vector<double>& mas) {
	for (auto& i : mas)
		out << i << " ";
	return out;
}

void write_iterations_information(double t0, const Result& result, const Function2D& u_true) {
	cout << "time: " << t0 << endl;
	cout << "residual: " << result.residual << endl;
	cout << "iterations: " << result.iterations << endl;
	cout << "answer:   " << result.answer.q << endl;

	lin_approx_t u_truly_approx = result.answer;
	for (int j = 0; j < u_truly_approx.x.size(); j++)
		u_truly_approx.q[j] = u_true(u_truly_approx.x[j], t0);
	cout << "shold be: " << u_truly_approx.q << endl;

	cout << "norm: " << norm(u_truly_approx.q, result.answer.q) << endl;

	cout << endl << "---------------------------" << endl;
}

void solveByTime(
	Function2D f,
	Function2D u_true,
	Function1D lambda,
	Function1D sigma,
	const lin_approx_t& u_start,
	const vector<double>& time_grid, 
	double eps,
	double maxiter
) {
	lin_approx_t u = u_start;

	for (int i = 0; i < time_grid.size()-1; ++i) {
		double t0 = time_grid[i];
		double t1 = time_grid[i+1];
		auto result = solveFixedPointIteration(f, u_true, lambda, sigma, u, t1-t0, t0, eps, maxiter);

		// Выводим мета-информацию
		write_iterations_information(t0, result, u_true);

		u = result.answer;
	}
}

int main() {
	auto u_true = [] (double x, double t) -> double { return 3*x*x + exp(t); };
	auto lambda = [] (double u) -> double { return u; };
	auto sigma = [] (double x) -> double { return 1; };
	auto f = calcRightPart(lambda, u_true, sigma);

	lin_approx_t u;
	make_grid(u.x, 0, 1, 10);
	for (int i = 0; i < u.x.size(); i++)
		u.q.push_back(u_true(u.x[i], 0));
		//u.q.push_back(u_true(u.x[i], 0)-0.01);
		//u.q.push_back(1);

	/*auto result = solveFixedPointIteration(f, u_true, lambda, sigma, u, 0.001, 0, 0.001, 1000);
	write_iterations_information(0, result, u_true); */

	vector<double> time;
	make_grid(time, 0, 1, 100);
	time.erase(time.begin());
	solveByTime(f, u_true, lambda, sigma, u, time, 1e-7, 100);

	system("pause");
}