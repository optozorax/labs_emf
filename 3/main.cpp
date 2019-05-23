#include "../2/lib.h"

struct Constants
{
	double sigma, lambda, xi, omega;
};

//-----------------------------------------------------------------------------
Function1D calcRightPartS(
	const Function1D& us,
	const Function1D& uc,
	const Constants& c
) {
	// fs = -lambda * div(grad us) - omega * sigma * uc - omega^2 * xi * uc
	return [=](double x) -> double {
		using namespace placeholders;
		auto divgrad = calcSecondDerivative(us);
		return -c.lambda * divgrad(x) - c.omega * c.sigma * uc(x) - c.omega * c.omega * c.xi * us(x);
	};
}

//-----------------------------------------------------------------------------
Function1D calcRightPartC(
	const Function1D& us,
	const Function1D& uc,
	const Constants& c
) {
	// fs = -lambda * div(grad uc) - omega * sigma * us - omega^2 * xi * uc
	return [=](double x) -> double {
		using namespace placeholders;
		auto divgrad = calcSecondDerivative(uc);
		return -c.lambda * divgrad(x) + c.omega * c.sigma * us(x) - c.omega * c.omega * c.xi * uc(x);
	};
}

//-----------------------------------------------------------------------------
double calcPIntegral(int i, int j, const lin_approx_t& u, const Constants& c) {
	auto f = [&] (double x) -> double {
		return c.lambda * u.basic_grad(x, i) * u.basic_grad(x, j) - c.omega * c.omega * c.xi * u.basic(x, i) * u.basic(x, j);
	};
	vector<double> X;
	make_grid(X, min(u.left(i), u.left(j)), max(u.right(i), u.right(j)), 20);
	return integral_gauss3(X, f);
}

//-----------------------------------------------------------------------------
double calcCIntegral(int i, int j, const lin_approx_t& u, const Constants& c) {
	auto f = [&] (double x) -> double {
		return u.basic(x, i) * u.basic(x, j);
	};
	vector<double> X;
	make_grid(X, min(u.left(i), u.left(j)), max(u.right(i), u.right(j)), 20);
	return c.omega * c.sigma * integral_gauss3(X, f);
}

//-----------------------------------------------------------------------------
Matrix calcLocalMatrix(int i, int j, const lin_approx_t& u, const Constants& c, bool isCalcLeftDown) {
	Matrix result(4, 4);
	double p11 = calcPIntegral(i, i, u, c);
	double c11 = calcCIntegral(i, i, u, c);
	double p12 = calcPIntegral(i, j, u, c);
	double c12 = calcCIntegral(i, j, u, c);
	double p21 = p12;
	double c21 = c12;
	double p22 = (isCalcLeftDown) ? calcPIntegral(j, j, u, c) : 1;
	double c22 = (isCalcLeftDown) ? calcCIntegral(j, j, u, c) : 1;

	result <<
		p11, -c11, p12, -c12,
		c11, p11,  c12, p12,
		p21, -c21, p22, -c22,
		c21, p21,  c22, p22;
	return result;
}

//-----------------------------------------------------------------------------
Matrix calcGlobalMatrix(const lin_approx_t& u, const Constants& c) {
	Matrix result(u.size() * 2, u.size() * 2);
	result.fill(0);

	for (int i = 0; i < u.size()-1; ++i) {
		auto l = calcLocalMatrix(i, i+1, u, c, i == u.size()-2);
		for (int x = 0; x < l.rows(); ++x) {
			for (int y = 0; y < l.cols(); ++y) {
				result(i*2 + x, i*2 + y) = l(x, y);
			}
		}
	}

	return result;
}

//-----------------------------------------------------------------------------
Vector calcB(
	const lin_approx_t& u, 
	const Constants& c, 
	const Function1D& fs,
	const Function1D& fc
) {
	Vector result(u.size() * 2);
	result.fill(0);
	for (int i = 0; i < u.size(); ++i) {
		auto funs = [&] (double x) -> double { return fs(x) * u.basic(x, i); };
		auto func = [&] (double x) -> double { return fc(x) * u.basic(x, i); };
		vector<double> X;
		make_grid(X, u.left(i), u.right(i), 10);
		result(i*2 + 0) = integral_gauss3(X, funs);
		result(i*2 + 1) = integral_gauss3(X, func);
	}
	return result;	
}

//-----------------------------------------------------------------------------
void setAnswer(
	lin_approx_t& us, 
	lin_approx_t& uc,
	const Vector& result 
) {
	us.q.resize(us.size());
	uc.q.resize(uc.size());
	for (int i = 0; i < us.size(); ++i) {
		us.q[i] = result(i*2 + 0);
		uc.q[i] = result(i*2 + 1);
	}
}

//-----------------------------------------------------------------------------
void setFirstBoundaryConditions(
	Matrix& A,
	Vector& b,
	const Function1D& us_true,
	const Function1D& uc_true,
	const lin_approx_t& u
) {
	auto clear_line = [&] (int line) {
		for (int i = 0; i < A.cols(); ++i) A(line, i) = 0;
	};

	clear_line(0); A(0, 0) = 1; b(0) = us_true(u.middle(0));
	clear_line(1); A(1, 1) = 1; b(1) = uc_true(u.middle(0));

	int end = A.cols()-1;
	clear_line(end-1); A(end-1, end-1) = 1; b(end-1) = us_true(u.middle(u.size()-1));
	clear_line(end); A(end, end) = 1; b(end) = uc_true(u.middle(u.size()-1));
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	cout.precision(3);

	Constants c;
	c.lambda = 32;
	c.omega = 100;
	c.xi = 10;
	c.sigma = 24;
	auto us_true = [] (double x) -> double { return 3*x*x + 2; };
	auto uc_true = [] (double x) -> double { return 6*x - x*x; };

	auto fs = calcRightPartS(us_true, uc_true, c);
	auto fc = calcRightPartC(us_true, uc_true, c);

	lin_approx_t us, uc;
	make_grid(us.x, 0, 1, 5);
	uc.x = us.x;

	auto A = calcGlobalMatrix(us, c);
	auto b = calcB(us, c, fs, fc);
	setFirstBoundaryConditions(A, b, us_true, uc_true, us);
	Eigen::JacobiSVD<Matrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Vector x = svd.solve(b);
	//cout << A << endl << b << endl << x << endl;

	setAnswer(us, uc, x);

	double residual = norm(us_true, us) + norm(uc_true, uc);

	cout << "answer s:    " << us.q << endl;
	cout << "should be s: " << calcTrulyApprox(us.x, us_true).q << endl;

	cout << "answer s:    " << uc.q << endl;
	cout << "should be s: " << calcTrulyApprox(uc.x, uc_true).q << endl;

	cout << "residual: " << residual << endl;

	system("pause");
}