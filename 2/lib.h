#pragma once

#include <functional>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

using namespace std;

//-----------------------------------------------------------------------------
// Типы данных
typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

typedef function<double(double)> Function1D;
typedef function<double(double, double)> Function2D;
typedef function<double(double, double, double)> Function3D;

//-----------------------------------------------------------------------------
// Интегралы
double integral_gauss3(const vector<double>& X, Function1D f);
void make_grid(vector<double>& X, double a, double b, long n, 
               Function1D move = [] (double x) -> double {return x;});

//-----------------------------------------------------------------------------
// Линейная базисная функция
double basicFunction(double left, double middle, double right, double x);
double basicFunctionGrad(double left, double middle, double right, double x);

//-----------------------------------------------------------------------------
// Производные
Function1D calcFirstDerivative(const Function1D& f);
Function1D calcSecondDerivative(const Function1D& f);

//-----------------------------------------------------------------------------
// Вывод информации на экран
string write_for_latex_double(double v, int precision);
ostream& operator<<(ostream& out, const vector<double>& mas);

//-----------------------------------------------------------------------------
// Линейная аппроксимация функции
class lin_approx_t
{
public:
	vector<double> q; /// Вектор весов
	vector<double> x; /// Массив положений каждого элемента

	int size(void) const;
	double left(int i) const;
	double middle(int i) const;
	double right(int i) const;

	double value(double pos) const; /// Получить значение функции аппроксимации в некоторой точке

	double basic(double pos, int i) const; /// Получить значение базовой функции под номером i в точке pos
	double basic_grad(double pos, int i) const ; /// Получить значение производной базовой функции под номером i в точке pos
};

//-----------------------------------------------------------------------------
// Результат работы итерационного процесса
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
ostream& operator<<(ostream& out, const Result::ExitType& exit);

//-----------------------------------------------------------------------------
// Сравнение результатов
lin_approx_t calcTrulyApprox(const vector<double>& grid, const Function2D& u_true, double time);
lin_approx_t calcTrulyApprox(const vector<double>& grid, const Function1D& u_true);
double norm(const vector<double>& a, const vector<double>& b); /// Считает норму между двумя векторами
double norm(const Function1D& f, const lin_approx_t& b); /// Считает норму между функцией и линейной аппроксимацией функции