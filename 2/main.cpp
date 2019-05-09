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

	double right(int i) const {
		if (i == x.size() - 1)
			return x.back()+0.00001;
		else
			return x[i+1];
	}

	double value(double pos) const {
		/*int start = std::distance(x.begin(), std::lower_bound(x.begin(), x.end(), pos))-1;
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
		}*/
		double sum = 0;
		//for (int i = start; i < std::min<int>(start+2, x.size()); ++i)
		for (int i = 0; i < x.size(); ++i)
			sum += q[i] * basic(pos, i);
		return sum;
	} /// Получить значение функции аппроксимации в некоторой точке

	double basic(double pos, int i) const {
		return basicFunction(left(i), x[i], right(i), pos);
	}

	double basic_grad(double pos, int i) const {
		return basicFunctionGrad(left(i), x[i], right(i), pos);	
	}
};

double calc_a1_integral(const lin_approx_t& lambda, const lin_approx_t& u, int i, int j) {
	if (std::abs(i - j) > 1)
		return 0;
	auto f = [&] (double x) -> double {
		return lambda.value(u.value(x)) * u.basic_grad(x, i) * u.basic_grad(x, j);
	};
	std::vector<double> X;
	make_grid(X, std::min(u.left(i), u.left(j))-1, std::max(u.right(i), u.right(j))+1, 100);
	return integral_gauss3(X, f);
}

double calc_a2_integral(const lin_approx_t& sigma, const lin_approx_t& u, int i, int j) {
	if (std::abs(i - j) > 1)
		return 0;
	auto f = [&] (double x) -> double {
		return sigma.value(x) * u.basic(x, i) * u.basic(x, j);
	};
	std::vector<double> X;
	make_grid(X, std::min(u.left(i), u.left(j)), std::max(u.right(i), u.right(j)), 5);
	return integral_gauss3(X, f);	
}

double calc_b1_integral(const lin_approx_t& f, const lin_approx_t& u, int i) {
	auto fun = [&] (double x) -> double {
		return f.value(x) * u.basic(x, i);
	};
	std::vector<double> X;
	make_grid(X, u.left(i)-1, u.right(i)+1, 100);
	return integral_gauss3(X, fun);	
}

double calc_b2_integral(const lin_approx_t& sigma, const lin_approx_t& u, const lin_approx_t& u_last_time, int i) {
	auto fun = [&] (double x) -> double {
		return sigma.value(x) * u_last_time.value(x) * u.basic(x, i);
	};
	std::vector<double> X;
	make_grid(X, u.left(i)-1, u.right(i)+1, 100);
	return integral_gauss3(X, fun);	
}

lin_approx_t approximate_function(Function1D lambda, const lin_approx_t& u) {
	lin_approx_t result = u;
	for (int i = 0; i < result.q.size(); ++i)
		result.q[i] = lambda(u.x[i]);
	return result;
}

Matrix calcA(const lin_approx_t& u_last, double dt, Function1D lambda, Function1D sigma) {
	auto lambda_approx = approximate_function(lambda, u_last);
	auto sigma_approx = approximate_function(sigma, u_last);

	Matrix result(u_last.size(), u_last.size());
	for (int i = 0; i < u_last.size(); ++i) {
		for (int j = 0; j < u_last.size(); ++j) {
			result(i, j) = calc_a1_integral(lambda_approx, u_last, i, j) + calc_a2_integral(sigma_approx, u_last, i, j) / dt;
		}
	}

	return result;
} /// Рассчитать матрицу A(q)
Vector calcB(const lin_approx_t& u_last, double dt, Function1D f, Function1D sigma, const lin_approx_t& u_last_time) {
	auto f_approx = approximate_function(f, u_last);
	auto sigma_approx = approximate_function(sigma, u_last);

	Vector result(u_last.size());
	for (int i = 0; i < u_last.size(); ++i) {
		result(i) = calc_b1_integral(f_approx, u_last, i) + calc_b2_integral(sigma_approx, u_last, u_last_time, i) / dt;
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
	int iterations;

	lin_approx_t answer;

	std::vector<Matrix> a_history;
	std::vector<Vector> b_history;
	std::vector<lin_approx_t> u_history;
};

Result solveFixedPointIteration(
	Function2D f,
	Function2D u_true,
	Function1D lambda,
	Function1D sigma,
	const lin_approx_t& u_last_time, 
	double dt,
	double time
) {
	lin_approx_t u;
	u.x = u_last_time.x;
	u.q.resize(u_last_time.x.size(), 1);

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
		result.a_history.push_back(A);
		result.b_history.push_back(b);
		result.u_history.push_back(u);

		if (result.iterations > 100) {
			result.exitType = Result::EXIT_ITERATIONS;
			break;
		}

		if ((A * x - b).norm()/b.norm() < 0.001) {
			result.exitType = Result::EXIT_RESIDUAL;
			break;
		}

		if ((x-last_x).norm() < 0.001) {
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
		const double h = 0.00001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}

Function2D calcRightPart(
	const Function1D& lambda, 
	const Function2D& u, 
	double sigma
) {
	// f = -div(lambda(u) * grad u) + sigma * du/dt
	return [=](double x, double t) -> double {
		using namespace std::placeholders;
		auto ut = calcFirstDerivative(std::bind(u, x, _1));
		auto ux = calcFirstDerivative(std::bind(u, _1, t));
		auto lambda_grad = [=](double x, double t) -> double {
			return lambda(u(x, t)) * ux(t);
		};
		auto div = calcFirstDerivative(std::bind(lambda_grad, _1, t));
		return -div(x) + sigma * ut(x);
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
		cout << i << " ";
	return out;
}


int main() {
	auto u_true = [] (double x, double t) -> double { return x; };
	auto lambda = [] (double u) -> double { return exp(exp(u)-1)+1; };
	auto sigma = [] (double x) -> double { return 0; };
	auto f = calcRightPart(lambda, u_true, 0);

	lin_approx_t u;
	make_grid(u.x, 0, 1, 10);
	for (int i = 0; i < u.x.size(); i++)
		u.q.push_back(u_true(u.x[i], 0));

	auto result = solveFixedPointIteration(f, u_true, lambda, sigma, u, 0.001, 0);

	cout << "iterations: " << result.iterations << endl;
	cout << "answer:   " << result.answer.q << endl;
	cout << "shold be: " << u.q << endl;
	cout << "norm: " << norm(u.q, result.answer.q) << endl;
	
	std::vector<double> grid;
	make_grid(grid, 0, 1, 100);
	cout << "integral norm: " << integral_gauss3(grid, [&](double x) -> double {
		return std::fabs(u.value(x)-result.answer.value(x));
	}) << endl;

	system("pause");
}