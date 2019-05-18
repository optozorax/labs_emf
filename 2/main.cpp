#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>

using namespace std;

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

typedef function<double(double)> Function1D;
typedef function<double(double, double)> Function2D;
typedef function<double(double, double, double)> Function3D;

//-----------------------------------------------------------------------------
double integral_gauss3(const vector<double>& X, Function1D f) {
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
void make_grid(vector<double>& X, double a, double b, long n) {
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

//-----------------------------------------------------------------------------
Function1D calcFirstDerivative(const Function1D& f) {
	return [f](double x) -> double {
		const double h = 0.001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}

//-----------------------------------------------------------------------------
Function2D calcRightPart(
	const Function1D& lambda, 
	const Function2D& u, 
	const Function1D& sigma
) {
	// f = -div(lambda(u) * grad u) + sigma * du/dt
	return [=](double x, double t) -> double {
		using namespace placeholders;
		auto ut = calcFirstDerivative(bind(u, x, _1));
		auto ux = calcFirstDerivative(bind(u, _1, t));
		auto lambda_grad = [lambda, ux, u](double x, double t) -> double {
			return lambda(u(x, t)) * ux(x);
		};
		auto div = calcFirstDerivative(bind(lambda_grad, _1, t));
		return -div(x) + sigma(x) * ut(x);
	};
}

//-----------------------------------------------------------------------------
double norm(const vector<double>& a, const vector<double>& b) {
	assert(a.size() == b.size());
	double sum = 0;
	for (int i = 0; i < a.size(); i++)
		sum += (a[i] - b[i])*(a[i] - b[i]);
	return sqrt(sum);
}

//-----------------------------------------------------------------------------
ostream& operator<<(ostream& out, const vector<double>& mas) {
	for (auto& i : mas)
		out << i << " ";
	return out;
}

//-----------------------------------------------------------------------------
string write_for_latex_double(double v, int precision) {
	int power = log(std::fabs(v)) / log(10.0);
	double value = v / pow(10.0, power);

	if (v != v) return "nan";

	if (v == 0) {
		power = 0;
		value = 0;
	}

	stringstream sout;
	sout.precision(precision);
	if (power == -1 || power == 0 || power == 1) {
		sout << v;
	} else {
		sout << value << "\\cdot 10^{" << power << "}";
	}

	return sout.str();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
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
		int start = distance(x.begin(), lower_bound(x.begin(), x.end(), pos))-1;
		if (start == -1) {
			if (fabs(pos-x[0]) < 0.00001)
				return q[0];
			else
				return 0;
		}
		if (start == x.size()) {
			if (fabs(pos-x.back()) < 0.00001)
				return q.back();
			else
				return 0;
		}
		double sum = 0;
		for (int i = start; i < min<int>(start+2, x.size()); ++i)
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

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
double calc_a1_integral(const Function1D& lambda, const lin_approx_t& u, int i, int j) {
	if (abs(i - j) > 1)
		return 0;
	auto f = [&] (double x) -> double {
		return lambda(u.value(x)) * u.basic_grad(x, i) * u.basic_grad(x, j);
	};
	vector<double> X;
	if (i != j)
		make_grid(X, min(u.middle(i), u.middle(j)), min(u.right(i), u.right(j)), 10);
	else
		make_grid(X, u.left(i), u.right(i), 10);
	return integral_gauss3(X, f);
}

//-----------------------------------------------------------------------------
double calc_a2_integral(const Function1D& sigma, const lin_approx_t& u, int i, int j) {
	if (abs(i - j) > 1)
		return 0;
	auto f = [&] (double x) -> double {
		return sigma(x) * u.basic(x, i) * u.basic(x, j);
	};
	vector<double> X;
	if (i != j)
		make_grid(X, min(u.middle(i), u.middle(j)), min(u.right(i), u.right(j)), 10);
	else
		make_grid(X, u.left(i), u.right(i), 10);
	return integral_gauss3(X, f);	
}

//-----------------------------------------------------------------------------
double calc_b1_integral(const Function1D& f, const lin_approx_t& u, int i) {
	auto fun = [&] (double x) -> double {
		return f(x) * u.basic(x, i);
	};
	vector<double> X;
	make_grid(X, u.left(i), u.right(i), 10);
	return integral_gauss3(X, fun);	
}

//-----------------------------------------------------------------------------
double calc_b2_integral(const Function1D& sigma, const lin_approx_t& u, const lin_approx_t& u_last_time, int i) {
	auto fun = [&] (double x) -> double {
		return sigma(x) * u_last_time.value(x) * u.basic(x, i);
	};
	vector<double> X;
	make_grid(X, u.left(i), u.right(i), 10);
	return integral_gauss3(X, fun);	
}

//-----------------------------------------------------------------------------
lin_approx_t approximate_function(Function1D lambda, const lin_approx_t& u) {
	lin_approx_t result = u;
	for (int i = 0; i < result.q.size(); ++i)
		result.q[i] = lambda(u.x[i]);
	return result;
}

//-----------------------------------------------------------------------------
Matrix calcA(const lin_approx_t& u_last, double dt, Function1D lambda, Function1D sigma) {
	Matrix result(u_last.size(), u_last.size());
	for (int i = 0; i < u_last.size()-1; ++i) {
		for (int j = 0; j < u_last.size(); ++j) {
			result(i, j) = calc_a1_integral(lambda, u_last, i, j) + calc_a2_integral(sigma, u_last, i, j) / dt;
		}
	}

	return result;
} /// Рассчитать матрицу A(q)

//-----------------------------------------------------------------------------
Vector calcB(const lin_approx_t& u_last, double dt, Function1D f, Function1D sigma, const lin_approx_t& u_last_time) {
	Vector result(u_last.size());
	for (int i = 0; i < u_last.size(); ++i) {
		result(i) = calc_b1_integral(f, u_last, i) + calc_b2_integral(sigma, u_last, u_last_time, i) / dt;
	}

	return result;
} /// Рассчитать вектор правой части b(q)

//-----------------------------------------------------------------------------
void write_first_boundary_conditions(Matrix& A, Vector& b, const lin_approx_t& u, Function2D u_true, double time) {
	A(0, 0) = 1;
	for (int i = 1; i < A.cols(); ++i) A(0, i) = 0;
	b(0) = u_true(u.x.front(), time);

	A(A.rows()-1, A.cols()-1) = 1;
	for (int i = 0; i < A.cols() - 1; ++i) A(A.rows()-1, i) = 0;
	b(b.rows()-1) = u_true(u.x.back(), time);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
ostream& operator<<(ostream& out, const Result::ExitType& exit) {
	switch (exit) {
		case Result::EXIT_RESIDUAL: out << "residual"; break;
		case Result::EXIT_STEP: out << "step"; break;
		case Result::EXIT_ITERATIONS: out << "iterations"; break;
		case Result::EXIT_ERROR: out << "error"; break;
	}
	return out;
}

//-----------------------------------------------------------------------------
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
		auto b = calcB(u, dt, bind(f, placeholders::_1, time), sigma, u_last_time);
		write_first_boundary_conditions(A, b, u, u_true, time);

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
			auto b1 = calcB(u, dt, bind(f, placeholders::_1, time), sigma, u_last_time);
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

lin_approx_t calcTrulyApprox(const vector<double>& grid, const Function2D& u_true, double time) {
	lin_approx_t u_truly_approx;
	u_truly_approx.x = grid;
	u_truly_approx.q = grid;
	for (int i = 0; i < u_truly_approx.x.size(); i++)
		u_truly_approx.q[i] = u_true(u_truly_approx.x[i], time);
	return u_truly_approx;
}

//-----------------------------------------------------------------------------
void write_iterations_information(double t0, const Result& result, const Function2D& u_true) {
	cout << "time: " << t0 << endl;
	cout << "residual: " << result.residual << endl;
	cout << "iterations: " << result.iterations << endl;
	cout << "answer:   " << result.answer.q << endl;

	lin_approx_t u_truly_approx = calcTrulyApprox(result.answer.x, u_true, t0);
	cout << "shold be: " << u_truly_approx.q << endl;

	cout << "norm: " << norm(u_truly_approx.q, result.answer.q) << endl;

	cout << endl << "---------------------------" << endl;
}

//-----------------------------------------------------------------------------
vector<Result> solveByTime(
	Function2D f,
	Function2D u_true,
	Function1D lambda,
	Function1D sigma,
	const lin_approx_t& u_start,
	const vector<double>& time_grid, 
	double eps,
	double maxiter
) {
	vector<Result> res;
	lin_approx_t u = u_start;

	for (int i = 0; i < time_grid.size()-1; ++i) {
		double t0 = time_grid[i];
		double t1 = time_grid[i+1];
		auto result = solveFixedPointIteration(f, u_true, lambda, sigma, u, t1-t0, t0, eps, maxiter);
		res.push_back(result);

		// Выводим мета-информацию
		// write_iterations_information(t0, result, u_true);

		u = result.answer;
	}

	return res;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<class A, class B, class C>
struct triple {
	A first;
	B second;
	C third;
};

//-----------------------------------------------------------------------------
void writeFirstInvestigation() {
	vector<triple<Function2D, string, bool>> us;
	us.push_back({[] (double x, double t) -> double { return 3*x + t; }, "$3x+t$", false});
	us.push_back({[] (double x, double t) -> double { return 2*x*x + t; }, "$2x^2+t$", false});
	us.push_back({[] (double x, double t) -> double { return x*x*x + t; }, "$x^3+t$", false});
	us.push_back({[] (double x, double t) -> double { return x*x*x*x + t; }, "$x^4+t$", false});
	us.push_back({[] (double x, double t) -> double { return exp(x) + t; }, "$e^x+t$", false});
	us.push_back({[] (double x, double t) -> double { return 3*x + t*t; }, "$3x+t^2$", true});
	us.push_back({[] (double x, double t) -> double { return 3*x + t*t*t; }, "$3x+t^3$", true});
	us.push_back({[] (double x, double t) -> double { return 3*x + exp(t); }, "$3x+e^t$", true});
	us.push_back({[] (double x, double t) -> double { return 3*x + sin(t); }, "$3x+sin(t)$", true});
	us.push_back({[] (double x, double t) -> double { return exp(x) + t*t; }, "$e^x+t^2$", true});
	us.push_back({[] (double x, double t) -> double { return exp(x) + t*t*t; }, "$e^x+t^3$", true});
	us.push_back({[] (double x, double t) -> double { return exp(x) + exp(t); }, "$e^x+e^t$", true});
	us.push_back({[] (double x, double t) -> double { return exp(x) + sin(t); }, "$e^x+sin(t)$", true});

	vector<pair<Function1D, string>> lambdas;
	lambdas.push_back({[] (double u) -> double { return 1; }, "$1$"});
	lambdas.push_back({[] (double u) -> double { return u; }, "$u$"});
	lambdas.push_back({[] (double u) -> double { return u*u; }, "$u^2$"});
	lambdas.push_back({[] (double u) -> double { return u*u + 1; }, "$u^2+1$"});
	lambdas.push_back({[] (double u) -> double { return u*u*u; }, "$u^3$"});
	lambdas.push_back({[] (double u) -> double { return u*u*u*u; }, "$u^4$"});
	lambdas.push_back({[] (double u) -> double { return exp(u); }, "$e^u$"});
	lambdas.push_back({[] (double u) -> double { return sin(u); }, "sinu"});

	double a = 1, b = 2, n = 10; // Характеристики сетки по пространству
	double at = 0, bt = 1, nt = 10; // Характеристики сетки по времени

	ofstream fout("first.txt");
	fout << "a\t";
	for (auto& i : lambdas) fout << i.second << ((i.second != lambdas.back().second) ? "\t" : "");
	fout << endl;

	vector<double> time;
	make_grid(time, at, bt, nt);
	time.erase(time.begin());
	time.push_back(time[1]-time[0] + bt);
	for (auto& i : us) {
		fout << i.second << "\t";
		auto u_true = i.first;
		auto sigma = [] (double x) -> double { return 1; };

		lin_approx_t u;
		make_grid(u.x, a, b, n);
		for (int j = 0; j < u.x.size(); j++)
			if (i.third)
				u.q.push_back(u_true(u.x[j], at));
			else
				u.q.push_back(1);

		lin_approx_t u_truly_approx = calcTrulyApprox(u.x, u_true, bt);

		for (auto& j : lambdas) {
			auto lambda = j.first;

			auto result = solveByTime(calcRightPart(lambda, u_true, sigma), u_true, lambda, sigma, u, time, 0.001, 500);

			int itersum = 0;
			for (auto& k : result) itersum += k.iterations;

			double residual = norm(u_truly_approx.q, result.back().answer.q);

			fout << "\\scalebox{.55}{\\tcell{$" << itersum << "$\\\\$" << write_for_latex_double(residual, 2) << "$}}" << ((j.second != lambdas.back().second) ? "\t" : "");
		}
		fout << endl;
	}

	fout.close();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	// auto u_true = [] (double x, double t) -> double { return 3*x*x + exp(t); };
	// auto lambda = [] (double u) -> double { return u; };
	// auto sigma = [] (double x) -> double { return 1; };
	// auto f = calcRightPart(lambda, u_true, sigma);

	// lin_approx_t u;
	// make_grid(u.x, 0, 1, 10);
	// for (int i = 0; i < u.x.size(); i++)
	// 	u.q.push_back(u_true(u.x[i], 0));
	// 	//u.q.push_back(u_true(u.x[i], 0)-0.01);
	// 	//u.q.push_back(1);

	// /*auto result = solveFixedPointIteration(f, u_true, lambda, sigma, u, 0.001, 0, 0.001, 1000);
	// write_iterations_information(0, result, u_true); */

	// vector<double> time;
	// make_grid(time, 0, 1, 100);
	// time.erase(time.begin());
	// solveByTime(f, u_true, lambda, sigma, u, time, 1e-7, 100);

	// system("pause");
	

	writeFirstInvestigation();

	//writeGridInvestigation();

	//system("pause");
}

/* Исследования:
	По вертикали - функция u, по горизонтали - функция лямбда, в клетках - число итераций и итоговая норма
		Функции: 
			3*x + t
			2*x*x + t
			x*x*x + t
			x*x*x*x + t
			exp(x) + t
			3*x + t*t
			3*x + t*t*t
			3*x + exp(t)
			3*x + sin(t)

		Лямбды:
			1
			u
			u*u
			u*u + 1
			u*u*u
			exp(u)

	Далее для какой-нибудь функции тестируется уменьшение сетки в два, четыре и т. д. раз

	Тестируется неравномерное разбиение по пространству для функций:
		x*x

		сетка:
		увеличивающаяся по коэффициенту, график для этого построить
		f^(-1)

		exp(x)
		сетка:
		коэфф
		1/exp(x)

		x*x*x

 */ 