#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <thread>
#include <future>
#include <chrono>

#include "lib.h"

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
		return -div(x) + sigma(x) * ut(t);
	};
}

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
EMatrix calcA(const lin_approx_t& u_last, double dt, Function1D lambda, Function1D sigma) {
	EMatrix result(u_last.size(), u_last.size());
	for (int i = 0; i < u_last.size()-1; ++i) {
		for (int j = 0; j < u_last.size(); ++j) {
			result(i, j) = calc_a1_integral(lambda, u_last, i, j) + calc_a2_integral(sigma, u_last, i, j) / dt;
		}
	}

	return result;
} /// Рассчитать матрицу A(q)

//-----------------------------------------------------------------------------
EVector calcB(const lin_approx_t& u_last, double dt, Function1D f, Function1D sigma, const lin_approx_t& u_last_time) {
	EVector result(u_last.size());
	for (int i = 0; i < u_last.size(); ++i) {
		result(i) = calc_b1_integral(f, u_last, i) + calc_b2_integral(sigma, u_last, u_last_time, i) / dt;
	}

	return result;
} /// Рассчитать вектор правой части b(q)

//-----------------------------------------------------------------------------
void write_first_boundary_conditions(EMatrix& A, EVector& b, const lin_approx_t& u, Function2D u_true, double time) {
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

	EVector last_x(u.q.size());
	for (int i = 0; i < u.q.size(); i++)
		last_x(i) = u.q[i];

	Result result;
	result.iterations = 0;
	while (true) {
		auto A = calcA(u, dt, lambda, sigma);
		auto b = calcB(u, dt, bind(f, placeholders::_1, time), sigma, u_last_time);
		write_first_boundary_conditions(A, b, u, u_true, time);

		Eigen::JacobiSVD<EMatrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		EVector x = svd.solve(b);
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

//-----------------------------------------------------------------------------
void write_iterations_information(double t0, const Result& result, const Function2D& u_true) {
	cout << "time: " << t0 << endl;
	cout << "residual: " << result.residual << endl;
	cout << "iterations: " << result.iterations << endl;
	cout << "answer:   " << result.answer.q << endl;

	//lin_approx_t u_truly_approx = calcTrulyApprox(result.answer.x, u_true, t0);
	//cout << "shold be: " << u_truly_approx.q << endl;

	cout << "norm: " << norm(bind(u_true, placeholders::_1, t0), result.answer) << endl;

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

	for (int i = 1; i < time_grid.size(); ++i) {
		double t0 = time_grid[i];
		double t1 = time_grid[i-1];
		auto result = solveFixedPointIteration(f, u_true, lambda, sigma, u, t0-t1, t0, eps, maxiter);
		res.push_back(result);

		// Выводим мета-информацию
		//write_iterations_information(t0, result, u_true);

		u = result.answer;
	}

	return res;
}

//-----------------------------------------------------------------------------
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
vector<triple<Function2D, string, bool>> us;
vector<pair<Function1D, string>> lambdas;
vector<Function2D> moves;
auto sigma = [] (double x) -> double { return 1; };

double a = 1, b = 2, n = 10; // Характеристики сетки по пространству
double at = 0, bt = 1, nt = 10; // Характеристики сетки по времени
double na = 0, nb = 1, nn = 1000; // Характеристики сетки по параметру t, в зависимости от которого меняются неравномерные сетки
int sa = 5, sb = 200; // Характеристики сетки по размеру

//-----------------------------------------------------------------------------
void init() {
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

	lambdas.push_back({[] (double u) -> double { return 1; }, "$1$"});
	lambdas.push_back({[] (double u) -> double { return u; }, "$u$"});
	lambdas.push_back({[] (double u) -> double { return u*u; }, "$u^2$"});
	lambdas.push_back({[] (double u) -> double { return u*u + 1; }, "$u^2+1$"});
	lambdas.push_back({[] (double u) -> double { return u*u*u; }, "$u^3$"});
	lambdas.push_back({[] (double u) -> double { return u*u*u*u; }, "$u^4$"});
	lambdas.push_back({[] (double u) -> double { return exp(u); }, "$e^u$"});
	lambdas.push_back({[] (double u) -> double { return sin(u); }, "sinu"});

	moves.push_back([] (double x, double t) -> double { return pow(x, t); });
	moves.push_back([] (double x, double t) -> double { return pow(x, 1.0/t); });
	moves.push_back([] (double x, double t) -> double { return (pow(t, x)-1.0)/(t-1.0); });
	moves.push_back([] (double x, double t) -> double { return (pow(1.0/t, x)-1.0)/(1.0/t-1.0); });
}

//-----------------------------------------------------------------------------
void writeFirstInvestigation() {
	ofstream fout("first.txt");
	fout << "a\t";
	for (auto& i : lambdas) fout << i.second << ((i.second != lambdas.back().second) ? "\t" : "");
	fout << endl;

	vector<double> time;
	make_grid(time, at, bt, nt);
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

		//lin_approx_t u_truly_approx = calcTrulyApprox(u.x, u_true, bt);

		for (auto& j : lambdas) {
			auto lambda = j.first;

			auto result = solveByTime(calcRightPart(lambda, u_true, sigma), u_true, lambda, sigma, u, time, 1e-7, 500);

			int itersum = 0;
			for (auto& k : result) itersum += k.iterations;

			double residual = norm(bind(u_true, placeholders::_1, bt), result.back().answer);

			fout << "\\scalebox{.55}{\\tcell{$" << itersum << "$\\\\$" << write_for_latex_double(residual, 2) << "$}}" << ((j.second != lambdas.back().second) ? "\t" : "");
		}
		fout << endl;
	}

	fout.close();
}

//-----------------------------------------------------------------------------
void writeGridInvestigation(
	const Function2D& u_true, string su,
	const string& file,
	const Function1D& lambda, string slambda) {
	// Делаем равномерную сетку по времени
	vector<double> time;
	make_grid(time, at, bt, nt);

	pair<int, double> resu = {-1, 0};

	auto test_grid = [&resu, time, u_true, lambda] (const vector<double>& grid) -> pair<int, double> {
		lin_approx_t u; u.x = grid;
		for (int j = 0; j < u.x.size(); j++) u.q.push_back(1);
		//lin_approx_t u_truly_approx = calcTrulyApprox(u.x, u_true, bt);

		auto result = solveByTime(calcRightPart(lambda, u_true, sigma), u_true, lambda, sigma, u, time, 0.001, 100);

		int itersum = 0; for (auto& k : result) itersum += k.iterations;
		double residual = norm(bind(u_true, placeholders::_1, bt), result.back().answer);
		if (residual > resu.second * 2 && resu.first != -1) residual = resu.second * 2;
		if (itersum > resu.first * 2 && resu.first != -1) itersum = resu.first * 2;
		return {itersum, residual};
	};

	ofstream fout(file);
	fout << "u = " << su << ", lambda = " << slambda << endl;

	// Получаем результаты для равномерной сетки
	vector<double> x_uniform;
	make_grid(x_uniform, a, b, n);
	resu = test_grid(x_uniform);

	fout << "t\tru\tiu";
	int counter = 0;
	for (auto& i : moves) {
		counter++;
		fout << "\tr" << counter << "\ti" << counter;
	}
	fout << endl;

	// Запускаем отложенные вычисления
	vector<std::future<pair<int, double>>> results;
	for (int i = 1; i <= nn; i++) {
		double t = na + (nb-na) * i/nn;

		for (auto& move : moves) {
			results.push_back(async([move, t, test_grid]()->pair<int, double>{
				vector<double> grid;
				make_grid(grid, a, b, n, bind(move, placeholders::_1, t));
				return test_grid(grid);
			}));
		}
	}

	auto it = results.begin();
	for (int i = 1; i <= nn; i++) {
		double t = na + (nb-na) * i/nn;
		cout << "\r" << 100 * t << "      ";
		fout << t << "\t" << resu.second << "\t" << resu.first;

		for (auto& move : moves) {
			auto value = it->get();
			fout << "\t" << value.second << "\t" << value.first;
			it++;
		}

		fout << endl;
	}
	cout << endl;

	fout.close();
}
//-----------------------------------------------------------------------------
void writeGridInvestigationTime(
	const Function2D& u_true, string su,
	const string& file,
	const Function1D& lambda, string slambda) {
	pair<int, double> resu = {-1, 0};

	auto test_grid = [u_true, lambda, resu] (const vector<double>& grid) -> pair<int, double> {
		lin_approx_t u; make_grid(u.x, a, b, n);
		for (int j = 0; j < u.x.size(); j++) u.q.push_back(u_true(u.x[j], at));
		//lin_approx_t u_truly_approx = calcTrulyApprox(u.x, u_true, bt);

		auto result = solveByTime(calcRightPart(lambda, u_true, sigma), u_true, lambda, sigma, u, grid, 0.001, 100);

		int itersum = 0; for (auto& k : result) itersum += k.iterations;
		double residual = norm(bind(u_true, placeholders::_1, bt), result.back().answer);
		if (residual > resu.second * 2 && resu.first != -1) residual = resu.second * 2;
		if (itersum > resu.first * 2 && resu.first != -1) itersum = resu.first * 2;
		return {itersum, residual};
	};

	ofstream fout(file);
	fout << "u = " << su << ", lambda = " << slambda << endl;

	// Получаем результаты для равномерной сетки
	vector<double> uniform;
	make_grid(uniform, at, bt, nt);
	resu = test_grid(uniform);

	fout << "t\tru\tiu";
	int counter = 0;
	for (auto& i : moves) {
		counter++;
		fout << "\tr" << counter << "\ti" << counter;
	}
	fout << endl;

	// Запускаем отложенные вычисления
	vector<std::future<pair<int, double>>> results;
	for (int i = 1; i <= nn; i++) {
		double t = na + (nb-na) * i/nn;

		for (auto& move : moves) {
			results.push_back(async([move, t, test_grid]()->pair<int, double>{
				vector<double> grid;
				make_grid(grid, at, bt, nt, bind(move, placeholders::_1, t));
				return test_grid(grid);
			}));
		}
	}

	auto it = results.begin();
	for (int i = 1; i <= nn; i++) {
		double t = na + (nb-na) * i/nn;
		cout << "\r" << 100 * t << "      ";
		fout << t << "\t" << resu.second << "\t" << resu.first;

		for (auto& move : moves) {
			auto value = it->get();
			fout << "\t" << value.second << "\t" << value.first;
			it++;
		}

		fout << endl;
	}
	cout << endl;

	for (int i = 1; i <= nn; i++) {
		double t = na + (nb-na) * i/nn;
		cout << "\r" << 100 * t << "      ";
		fout << t << "\t" << resu.second << "\t" << resu.first;

		vector<double> grid;
		for (auto& move : moves) {
			make_grid(grid, at, bt, nt, bind(move, placeholders::_1, t));
			auto res = test_grid(grid);
			fout << "\t" << res.second << "\t" << res.first;
		}

		fout << endl;
	}

	cout << endl;

	fout.close();
}

//-----------------------------------------------------------------------------
void writeGridInvestigationSize(
	const Function2D& u_true, string su,
	const string& file,
	const Function1D& lambda, string slambda) {
	ofstream fout(file);
	fout << "u = " << su << ", lambda = " << slambda << endl;

	// Делаем равномерную сетку по времени
	vector<double> time;
	make_grid(time, at, bt, nt);

	fout << "size\titerations\tresidual" << endl;
	for (int size = sa; size < sb; ++size) {
		lin_approx_t u;
		make_grid(u.x, a, b, size);
		for (int j = 0; j < u.x.size(); j++) u.q.push_back(u_true(u.x[j], at));
		//lin_approx_t u_truly_approx = calcTrulyApprox(u.x, u_true, bt);

		auto result = solveByTime(calcRightPart(lambda, u_true, sigma), u_true, lambda, sigma, u, time, 0.001, 100);

		int itersum = 0; for (auto& k : result) itersum += k.iterations;
		double residual = norm(bind(u_true, placeholders::_1, bt), result.back().answer);

		residual /= size; // Важно! Так как если не разделить, то расстояние между каждым узлом будет в итоге уменьшаться, но раз мы их делаем всё больше и больше, то в сумме оно будет увеличиваться, и в итоге точность будет расти.

		fout << size << "\t" << itersum << "\t" << residual << endl;
	}

	fout.close();
}

//-----------------------------------------------------------------------------
void writeGridInvestigationSizeTime(
	const Function2D& u_true, string su,
	const string& file,
	const Function1D& lambda, string slambda) {
	ofstream fout(file);
	fout << "u = " << su << ", lambda = " << slambda << endl;

	lin_approx_t u;
	make_grid(u.x, a, b, n);
	for (int j = 0; j < u.x.size(); j++) u.q.push_back(u_true(u.x[j], at));
	//lin_approx_t u_truly_approx = calcTrulyApprox(u.x, u_true, bt);

	fout << "size\titerations\tresidual" << endl;
	for (int size = sa; size < sb; ++size) {
		vector<double> time;
		make_grid(time, at, bt, size);

		auto result = solveByTime(calcRightPart(lambda, u_true, sigma), u_true, lambda, sigma, u, time, 0.001, 100);

		int itersum = 0; for (auto& k : result) itersum += k.iterations;
		double residual = norm(bind(u_true, placeholders::_1, bt), result.back().answer);

		residual /= size; // Важно! Так как если не разделить, то расстояние между каждым узлом будет в итоге уменьшаться, но раз мы их делаем всё больше и больше, то в сумме оно будет увеличиваться, и в итоге точность будет расти.

		fout << size << "\t" << itersum << "\t" << residual << endl;
	}

	fout.close();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	init();

	vector<thread> t;

	/*t.emplace_back([](){
	writeFirstInvestigation();
	});*/

	t.emplace_back([](){
	writeGridInvestigation(
		[] (double x, double t) -> double { return x*x*x*x + t; }, "$x^4+t$", 
		"x4_space.txt",
		[] (double u) -> double { return u*u; }, "$u^2$"
	);
	});

	t.emplace_back([](){
	writeGridInvestigation(
		[] (double x, double t) -> double { return exp(x) + t; }, "$exp(x)+t$", 
		"expx_space.txt",
		[] (double u) -> double { return u*u; }, "$u^2$"
	);
	});

	t.emplace_back([](){
	writeGridInvestigationTime(
		[] (double x, double t) -> double { return exp(x) + t*t*t; }, "$e^x+t^3$", 
		"t3_time.txt",
		[] (double u) -> double { return u; }, "$u$"
	);
	});

	t.emplace_back([](){
	writeGridInvestigationTime(
		[] (double x, double t) -> double { return exp(x) + exp(t); }, "$e^x+e^t$", 
		"expt_time.txt",
		[] (double u) -> double { return u; }, "$u$"
	);
	});

	/*t.emplace_back([](){
	writeGridInvestigationSize(
		[] (double x, double t) -> double { return exp(x) + exp(t); }, "$e^x+e^t$", 
		"expx_expt_size_space.txt",
		[] (double u) -> double { return u; }, "$u$"
	);
	});

	t.emplace_back([](){
	writeGridInvestigationSizeTime(
		[] (double x, double t) -> double { return exp(x) + exp(t); }, "$e^x+e^t$", 
		"expx_expt_size_time.txt",
		[] (double u) -> double { return u; }, "$u$"
	);
	});*/

	auto start = chrono::high_resolution_clock::now();
	for (auto& i : t) i.join();
	auto end = chrono::high_resolution_clock::now();
	double time = chrono::duration_cast<chrono::seconds>(end-start).count();
	cout << "time: " << time << endl;
	system("pause");
}