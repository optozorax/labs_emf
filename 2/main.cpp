#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <Eigen/Dense>

using namespace std;

typedef Eigen::MatrixXd Matrix;

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
		return lambda.value(x) * u.basic_grad(x, i) * u.basic_grad(x, j);
	};
	std::vector<double> X;
	make_grid(X, std::min(u.left(i), u.left(j))-1, std::max(u.right(i), u.right(j))+1, 100);
	return integral_gauss3(X, f);
}

double calc_a2_integral(const lin_approx_t& sigma, const lin_approx_t& u, int i, int j) {
	//if (std::abs(i - j) > 1)
		//return 0;
	auto f = [&] (double x) -> double {
		return sigma.value(x) * u.basic(x, i) * u.basic(x, j);
	};
	std::vector<double> X;
	make_grid(X, std::min(u.left(i), u.left(j)), std::max(u.right(i), u.right(j)), 5);
	return integral_gauss3(X, f);	
}

double calc_b_integral(const lin_approx_t& f, const lin_approx_t& u, int i) {
	//if (std::abs(i - j) > 1)
	//return 0;
	auto fun = [&] (double x) -> double {
		return f.value(x) * u.basic(x, i);
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
	return {};
} /// Рассчитать матрицу A(q)
//Vector calcB(const Vector& q); /// Рассчитать вектор правой части b(q)

struct Result
{
	enum ExitType {
		EXIT_RESIDUAL,
		EXIT_STEP,
		EXIT_ITERATIONS,
		EXIT_ERROR
	} exitType;
	int iterations;

	
};

/*Result solveFixedPointIteration(
	const Vector& q0, 
	double eps
);*/ /// Решить уравнение методом простой итерации

double lambda(double ut) { return 1; } 

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

int main() {
	lin_approx_t u;
	make_grid(u.x, 0, 1, 10);
	for (int i = 0; i < u.x.size(); i++)
		u.q.push_back(i/10.0);

	auto l = [] (double x) -> double {return x*x + 1;};
	auto u_true = [] (double x, double t) -> double {return x;};

	auto lambda = approximate_function(l, u);

	for (int i = 0; i < u.x.size(); ++i) {
		for (int j = 0; j < u.x.size(); ++j) {
			cout << std::setw(10) << calc_a1_integral(lambda, u, i, j);
		}
		cout << endl;
	}

	auto f_true = std::bind(calcRightPart(l, u_true, 0), std::placeholders::_1, 0);
	auto f = approximate_function(f_true, u);

	for (int i = 0; i < u.x.size(); ++i) {
		cout << calc_b_integral(f, u, i) << endl;
	}

	cout << fixed;
	for (double i = 0; i < 1; i+=0.01)
		cout << "(" << i << ", " << f.value(i) * u.basic(i, 3) << ")" << endl;

	system("pause");
}