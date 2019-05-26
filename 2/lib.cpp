#include "lib.h"

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
void make_grid(vector<double>& X, double a, double b, long n, Function1D move) {
	X.clear();
	double size = b-a;
	for (double i = 0; i <= n; ++i)
		X.push_back(a + move(i/double(n)) * size);
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
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
Function1D calcFirstDerivative(const Function1D& f) {
	return [f](double x) -> double {
		const double h = 0.001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}

//-----------------------------------------------------------------------------
Function1D calcSecondDerivative(const Function1D& f) {
	return [f](double x) -> double {
		const double h = 0.0001;
		return (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h))/(12*h*h);
	};
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

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
int lin_approx_t::size(void) const {
	return x.size();
}

//-----------------------------------------------------------------------------
double lin_approx_t::left(int i) const {
	if (i == 0)
		return x[0]-0.00001;
	else
		return x[i-1];
}

//-----------------------------------------------------------------------------
double lin_approx_t::middle(int i) const {
	return x[i];
}

//-----------------------------------------------------------------------------
double lin_approx_t::right(int i) const {
	if (i == x.size() - 1)
		return x.back()+0.00001;
	else
		return x[i+1];
}

//-----------------------------------------------------------------------------
double lin_approx_t::value(double pos) const {
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
		sum += q[i] * basic(pos, i);
	return sum;
}

//-----------------------------------------------------------------------------
double lin_approx_t::basic(double pos, int i) const {
	return basicFunction(left(i), middle(i), right(i), pos);
}

//-----------------------------------------------------------------------------
double lin_approx_t::basic_grad(double pos, int i) const {
	return basicFunctionGrad(left(i), middle(i), right(i), pos);	
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

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
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
lin_approx_t calcTrulyApprox(const vector<double>& grid, const Function2D& u_true, double time) {
	lin_approx_t u_truly_approx;
	u_truly_approx.x = grid;
	u_truly_approx.q = grid;
	for (int i = 0; i < u_truly_approx.x.size(); i++)
		u_truly_approx.q[i] = u_true(u_truly_approx.x[i], time);
	return u_truly_approx;
}

//-----------------------------------------------------------------------------
lin_approx_t calcTrulyApprox(const vector<double>& grid, const Function1D& u_true) {
	lin_approx_t u_truly_approx;
	u_truly_approx.x = grid;
	u_truly_approx.q = grid;
	for (int i = 0; i < u_truly_approx.x.size(); i++)
		u_truly_approx.q[i] = u_true(u_truly_approx.x[i]);
	return u_truly_approx;
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
double norm(const Function1D& f, const lin_approx_t& b) {
	vector<double> grid;
	make_grid(grid, b.middle(0), b.middle(b.x.size()-1), b.size() * 10);
	return integral_gauss3(grid, [f, &b] (double x) -> double {
		return fabs(f(x) - b.value(x));
	});
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
Function1D getMove0(double t, double coef) {
	return getMove1(1.0/t, coef);
}

//-----------------------------------------------------------------------------
Function1D getMove1(double t, double coef) {
	t = pow(t, coef);
	return [t] (double x) -> double {
		return (1.0-pow(t, x))/(1.0-t);
	};
}