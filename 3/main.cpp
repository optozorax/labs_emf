#include <iomanip>

#include "../2/lib.h"
#include <diagonal.h>

//-----------------------------------------------------------------------------
double operator*(const vector<double>& a, const vector<double>& b) {
	double sum = 0;
	for (int i = 0; i < a.size(); ++i)
		sum += a[i] * b[i];
	return sum;
}

//-----------------------------------------------------------------------------
double length(const vector<double>& mas) {
	return sqrt(mas*mas);
}

//-----------------------------------------------------------------------------
vector<double> to(const Vector& a) {
	vector<double> result(a.size());
	for (int i = 0; i < a.size(); ++i)
		result[i] = a(i);
	return result;
}

//-----------------------------------------------------------------------------
Vector to(const vector<double>& a) {
	Vector result(a.size());
	for (int i = 0; i < a.size(); ++i)
		result(i) = a[i];
	return result;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
struct matrix
{
	int n;
	vector<double> d, l, u;
	vector<int> i, j;

	void init(int n1) {
		n = n1;
		d.clear();
		l.clear();
		u.clear();
		i.clear();
		j.clear();
		d.resize(n);
		i.resize(n+1, 0);
	}

	void toDense(Matrix& m) const {
		m.resize(n, n, 0);
		for (int i = 0; i < n; ++i) {
			m(i, i) = d[i];
			for (int j = 0; j < lineElemCount(i); ++j) {
				m(i, lineElemRow(i, j)) = l[lineElemStart(i) + j];
				m(lineElemRow(i, j), i) = u[lineElemStart(i) + j];
			}
		}
	}

	int lineElemStart(int line) const { 
		return i[line]; 
	}
	int lineStart(int line) const { 
		return j[lineElemStart(line)]; 
	}
	int lineSize(int line) const { 
		return line - lineStart(line); 
	} 
	int lineElemRow(int line, int elem) const { 
		return j[lineElemStart(line) + elem]; 
	}
	int lineElemCount(int line) const { 
		return i[line+1]-i[line]; 
	}
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

void lu_decompose(const matrix& a, matrix& lu) {
	lu = a;
	for (int i = 0; i < lu.n; ++i) {
		// Заполняем нижний треугольник
		int line_start = lu.lineElemStart(i);
		int line_end = lu.lineElemStart(i+1);
		for (int j = line_start; j < line_end; ++j) {
			double sum = 0;

			int row = lu.j[j];
			int row_start = lu.lineElemStart(row);
			int row_end = lu.lineElemStart(row+1);

			int kl = line_start;
			int ku = row_start;
			
			while (kl < j && ku < row_end) {
				if (lu.j[kl] == lu.j[ku]) { // Совпадают столбцы
					sum += lu.l[kl] * lu.u[ku];
					ku++;
					kl++;
				} else if (lu.j[kl] < lu.j[ku]) {
					kl++;
				} else {
					ku++;
				}
			}

			lu.l[j] = (a.l[j] - sum) / lu.d[row];
		}

		// Заполняем верхний треугольник
		int row_start = lu.lineElemStart(i);
		int row_end = lu.lineElemStart(i+1);
		for (int j = line_start; j < line_end; ++j) {
			double sum = 0;
			
			int line = lu.j[j];
			int line_start = lu.lineElemStart(line);
			int line_end = lu.lineElemStart(line+1);

			int kl = line_start;
			int ku = row_start;
			
			while (kl < line_end && ku < j) {
				if (lu.j[kl] == lu.j[ku]) { // Совпадают столбцы
					sum += lu.l[kl] * lu.u[ku];
					ku++;
					kl++;
				} else if (lu.j[kl] < lu.j[ku]) {
					kl++;
				} else {
					ku++;
				}
			}

			lu.u[j] = (a.u[j] - sum) / lu.d[line];
		}

		// Расчитываем диагональный элемент
		double sum = 0;
		int line_row_start = lu.lineElemStart(i);
		int line_row_end = lu.lineElemStart(i+1);
		for (int j = line_row_start; j < line_row_end; ++j)
			sum += lu.l[j] * lu.u[j];

		lu.d[i] = sqrt(a.d[i] - sum);
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void mul(const matrix& a, vector<double>& x_y) {
	vector<double> result(a.n, 0);

	for (int i = 0; i < a.n; ++i) {
		int start = a.lineElemStart(i);
		int size = a.lineElemCount(i);
		for (int j = 0; j < size; j++) {
			result[i] += a.l[start + j] * x_y[a.lineElemRow(i, j)];
			result[a.lineElemRow(i, j)] += a.u[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < a.n; ++i)
		result[i] += a.d[i] * x_y[i];

	x_y = result;
}

//-----------------------------------------------------------------------------
void mul_t(const matrix& a, vector<double>& x_y) {
	vector<double> result(a.n, 0);

	for (int i = 0; i < a.n; ++i) {
		int start = a.lineElemStart(i);
		int size = a.lineElemCount(i);
		for (int j = 0; j < size; j++) {
			result[i] += a.u[start + j] * x_y[a.lineElemRow(i, j)];
			result[a.lineElemRow(i, j)] += a.l[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < a.n; ++i)
		result[i] += a.d[i] * x_y[i];

	x_y = result;
}

//-----------------------------------------------------------------------------
void mul_l_invert_t(const matrix& l, vector<double>& y_x) {
	for (int i = l.n - 1; i >= 0; i--) {
		int start = l.lineElemStart(i);
		int size = l.lineElemCount(i);

		y_x[i] /= l.d[i];
		for (int j = 0; j < size; ++j)
			y_x[l.lineElemRow(i, j)] -= y_x[i] * l.l[start + j];
	}
}

//-----------------------------------------------------------------------------
void mul_u_invert_t(const matrix& u, vector<double>& y_x) {
	for (int i = 0; i < u.n; ++i) {
		int start = u.lineElemStart(i);
		int size = u.lineElemCount(i);

		sumreal sum = 0;
		for (int j = 0; j < size; ++j)
			sum += u.u[start + j] * y_x[u.lineElemRow(i, j)];
		y_x[i] = (y_x[i] - sum) / u.d[i];
	}
}

//-----------------------------------------------------------------------------
void mul_l_invert(const matrix& l, vector<double>& y_x) {
	for (int i = 0; i < l.n; ++i) {
		int start = l.lineElemStart(i);
		int size = l.lineElemCount(i);

		sumreal sum = 0;
		for (int j = 0; j < size; ++j)
			sum += l.l[start + j] * y_x[l.lineElemRow(i, j)];
		y_x[i] = (y_x[i] - sum) / l.d[i];
	}
}

//-----------------------------------------------------------------------------
void mul_u_invert(const matrix& u, vector<double>& y_x) {
	for (int i = u.n-1; i >= 0; i--) {
		int start = u.lineElemStart(i);
		int size = u.lineElemCount(i);

		y_x[i] /= u.d[i];
		for (int j = 0; j < size; ++j)
			y_x[u.lineElemRow(i, j)] -= y_x[i] * u.u[start + j];
	}
}

//-----------------------------------------------------------------------------
void mul_u(const matrix& u, vector<double>& x_y) {
	vector<double> result(u.n, 0);

	for (int i = 0; i < u.n; ++i) {
		int start = u.lineElemStart(i);
		int size = u.lineElemCount(i);
		for (int j = 0; j < size; j++) {
			result[u.lineElemRow(i, j)] += u.u[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < u.n; ++i)
		result[i] += u.d[i] * x_y[i];

	x_y = result;
}

//-----------------------------------------------------------------------------
void mul(const vector<double>& d, vector<double>& x_y) {
	for (int i = 0; i < d.size(); i++)
		x_y[i] *= d[i];
}

//-----------------------------------------------------------------------------
void mul_invert(const vector<double>& d, vector<double>& x_y) {
	for (int i = 0; i < d.size(); i++)
		x_y[i] /= d[i];
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
class SLAU
{
public:

//-----------------------------------------------------------------------------
pair<int, double> los2() {
	lu_decompose(a, lu);
	x.clear();
	x.resize(n, 0);

	r = x;
	mul(a, r);
	for (int i = 0; i < n; i++)
		r[i] = f[i] - r[i];
	mul_l_invert(lu, r);

	z = r;
	mul_u_invert(lu, z);

	p = z;
	mul(a, p);
	mul_l_invert(lu, p);

	double flen = sqrt(f*f);
	double residual;

	int i = 0;
	while (true) {
		double pp = p*p;
		double alpha = (p*r) / pp;
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * z[i];
			r[i] -= alpha * p[i];
		}
		t1 = r;
		mul_u_invert(lu, t1);
		t2 = t1;
		mul(a, t2);
		mul_l_invert(lu, t2);
		double beta = -(p*t2) / pp;
		for (int i = 0; i < n; ++i) {
			z[i] = t1[i] + beta * z[i];
			p[i] = t2[i] + beta * p[i];
		}
		residual = length(r) / flen;
		i++;

		if (is_log) cout << "Iteration: " << setw(4) << i << ", Residual: " << setw(20) << setprecision(16) << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}

	return {i, residual};
}

//-----------------------------------------------------------------------------
pair<int, double> bsg_stab_lu() {
	lu_decompose(a, lu);
	x.clear();
	x.resize(n, 0);
	vector<double> r0(n, 0);
	vector<double> y = x;
	mul(a, y);

	r = x;
	mul(a, r);
	for (int i = 0; i < n; i++)
		r[i] = f[i] - r[i];
	mul_l_invert(lu, r);

	z = r;
	mul_u_invert(lu, z);

	r0 = r; // r0 - это r0
	p = r;

	double flen = sqrt(f*f);
	double residual;

	int i = 0;
	while (true) {
		t1 = z;
		mul_u_invert(lu, t1);
		mul(a, t1);
		mul_l_invert(lu, t1);
 		double rr0 = r*r0;
		double alpha = (rr0) / (t1*r0);
		for (int i = 0; i < n; ++i)
			p[i] = r[i] - alpha * t1[i];

		t2 = p;
		mul_u_invert(lu, t2);
		mul(a, t2);
		mul_l_invert(lu, t2);
		double gamma = (p*t2) / (t2*t2);
		for (int i = 0; i < n; ++i) {
			y[i] = y[i] + alpha * z[i] + gamma * p[i];
			r[i] = p[i] - gamma * t2[i];
		}

		double beta = alpha*(r*r0)/(gamma * rr0);
		for (int i = 0; i < n; ++i)
			z[i] = r[i] + beta * z[i] - beta * gamma * t1[i];
		x = y;
		mul_u_invert(lu, x);

		residual = length(r) / flen;
		i++;

		if (is_log) cout << "Iteration: " << setw(4) << i << ", Residual: " << setw(20) << setprecision(16) << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}

	return {i, residual};
}

int n, maxiter;
double eps;
matrix a, lu;
vector<double> f;
vector<double> r, z, p;
vector<double> x, t1, t2;
bool is_log;

};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

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
	// fs = -lambda * div(grad uc) + omega * sigma * us - omega^2 * xi * uc
	return [=](double x) -> double {
		using namespace placeholders;
		auto divgrad = calcSecondDerivative(uc);
		return -c.lambda * divgrad(x) + c.omega * c.sigma * us(x) - c.omega * c.omega * c.xi * uc(x);
	};
}

const int count_integral = 50;

//-----------------------------------------------------------------------------
double calcPIntegral(int i, int j, const lin_approx_t& u, const Constants& c) {
	auto f = [&] (double x) -> double {
		return c.lambda * u.basic_grad(x, i) * u.basic_grad(x, j) - c.omega * c.omega * c.xi * u.basic(x, i) * u.basic(x, j);
	};
	vector<double> X;
	if (i == j)
		make_grid(X, u.left(i), u.right(i), count_integral);
	else
		make_grid(X, min(u.middle(i), u.middle(j)), min(u.right(i), u.right(j)), count_integral);
	return integral_gauss3(X, f);
}

//-----------------------------------------------------------------------------
double calcCIntegral(int i, int j, const lin_approx_t& u, const Constants& c) {
	auto f = [&] (double x) -> double {
		return u.basic(x, i) * u.basic(x, j);
	};
	vector<double> X;
	if (i == j)
		make_grid(X, u.left(i), u.right(i), count_integral);
	else
		make_grid(X, min(u.middle(i), u.middle(j)), min(u.right(i), u.right(j)), count_integral);
	return c.omega * c.sigma * integral_gauss3(X, f);
}

//-----------------------------------------------------------------------------
EMatrix calcLocalMatrix(int i, int j, const lin_approx_t& u, const Constants& c, bool isCalcLeftDown) {
	EMatrix result(4, 4);
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
EMatrix calcGlobalMatrix(const lin_approx_t& u, const Constants& c) {
	EMatrix result(u.size() * 2, u.size() * 2);
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
matrix calcGlobalMatrixProfile(const lin_approx_t& u, const Constants& c) {
	matrix result;
	result.n = u.size() * 2;
	result.d.resize(result.n, 0);
	result.i.resize(result.n+1, -1);

	int counter = 0;
	for (int i = 0; i < u.size()-1; ++i) {
		auto l = calcLocalMatrix(i, i+1, u, c, true);
		for (int y = 0; y < l.cols(); ++y) {
			if (result.i[i*2 + y] == -1) {
				result.i[i*2 + y] = counter;
			}
			if ((i != 0 && y > 1) || (i == 0)) {
				for (int x = 0; x <= y; ++x) {
					if (x == y) {
						result.d[i*2 + x] = l(x, y);
					} else {
						result.j.push_back(i*2 + x);
						result.u.push_back(l(x, y));
						result.l.push_back(l(y, x));
						counter++;
					}
				}
			}
		}
	}
	result.i.back() = counter;

	return result;
}

//-----------------------------------------------------------------------------
template<class V>
V calcB(
	const lin_approx_t& u, 
	const Constants& c, 
	const Function1D& fs,
	const Function1D& fc
) {
	V result(u.size() * 2);
	for (int i = 0; i < u.size(); ++i) {
		auto funs = [&] (double x) -> double { return fs(x) * u.basic(x, i); };
		auto func = [&] (double x) -> double { return fc(x) * u.basic(x, i); };
		vector<double> X;
		make_grid(X, u.left(i), u.right(i), count_integral);
		result(i*2 + 0) = integral_gauss3(X, funs);
		result(i*2 + 1) = integral_gauss3(X, func);
	}
	return result;	
}

//-----------------------------------------------------------------------------
template<class V>
void setAnswer(
	lin_approx_t& us, 
	lin_approx_t& uc,
	const V& result 
) {
	us.q.resize(us.size());
	uc.q.resize(uc.size());
	for (int i = 0; i < us.size(); ++i) {
		us.q[i] = result(i*2 + 0);
		uc.q[i] = result(i*2 + 1);
	}
}

//-----------------------------------------------------------------------------
MatrixDiagonal calcGlobalMatrixDiag(const lin_approx_t& u, const Constants& c) {
	Diagonal d(u.size() * 2);
	vector<int> format = {0, -3, -2, -1, 1, 2, 3};
	MatrixDiagonal result(u.size() * 2, format);
	for (int i = 0; i < u.size()-1; ++i) {
		auto l = calcLocalMatrix(i, i+1, u, c, i == u.size()-2);
		for (int x = 0; x < l.rows(); ++x) {
			for (int y = 0; y < l.cols(); ++y) {
				int line = i*2 + x;
				int row = i*2 + y;
				int diag = d.calcDiag_byLR(line, row);
				int pos = d.calcPos_byLR(line, row);
				diag = distance(format.begin(), find(format.begin(), format.end(), diag));
				result.begin(diag)[pos] = l(x, y);
			}
		}
	}

	return result;
}

//-----------------------------------------------------------------------------
void setFirstBoundaryConditions(
	EMatrix& A,
	EVector& b,
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
	clear_line(end);   A(end, end) = 1;     b(end) = uc_true(u.middle(u.size()-1));
}

//-----------------------------------------------------------------------------
void setFirstBoundaryConditions(
	MatrixDiagonal& A,
	Vector& b,
	const Function1D& us_true,
	const Function1D& uc_true,
	const lin_approx_t& u
) {
	auto clear_line = [&] (int line) {
		matrix_diagonal_line_iterator it(A.dimension(), A.getFormat(), false);
		for (; !it.isEnd(); ++it) 
			for (; !it.isLineEnd(); ++it)
				if (it.i == line)
					A.begin(it.dn)[it.di] = 0;
	};

	clear_line(0); A.begin(0)[0] = 1; b(0) = us_true(u.middle(0));
	clear_line(1); A.begin(0)[1] = 1; b(1) = uc_true(u.middle(0));

	int end = A.dimension()-1;
	clear_line(end-1); A.begin(0)[end-1] = 1; b(end-1) = us_true(u.middle(u.size()-1));
	clear_line(end);   A.begin(0)[end] = 1;   b(end) = uc_true(u.middle(u.size()-1));
}

//-----------------------------------------------------------------------------
void setFirstBoundaryConditions(
	matrix& A,
	vector<double>& b,
	const Function1D& us_true,
	const Function1D& uc_true,
	const lin_approx_t& u
) {
	auto clear_lines = [&] (vector<int> lines) {
		for (auto& i : lines) A.d[i] = 0;

		for (int i = 0; i < A.i.size()-1; i++) {
			for (int pj = A.i[i]; pj < A.i[i+1]; pj++) {
				int j = A.j[pj];
				if (find(lines.begin(), lines.end(), j) != lines.end()) A.u[pj] = 0;
				if (find(lines.begin(), lines.end(), i) != lines.end()) A.l[pj] = 0;
			}
		}
	};

	int end = A.n-1;
	clear_lines({0, 1, end-1, end});

	A.d[0] = 1; b[0] = us_true(u.middle(0));
	A.d[1] = 1; b[1] = uc_true(u.middle(0));

	A.d[end-1] = 1; b[end-1] = us_true(u.middle(u.size()-1));
	A.d[end] = 1;   b[end] = uc_true(u.middle(u.size()-1));
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main() {
	cout.precision(3);

	Constants c;
	/*c.lambda = 32;
	c.omega = 100;
	c.xi = 10;
	c.sigma = 24;*/
	c.lambda = 1;
	c.omega = 1;
	c.xi = 1;
	c.sigma = 1;
	auto us_true = [] (double x) -> double { return 3*x*x*x*x + 2*exp(x); };
	auto uc_true = [] (double x) -> double { return 6*x - pow(x, exp(x)); };

	auto fs = calcRightPartS(us_true, uc_true, c);
	auto fc = calcRightPartC(us_true, uc_true, c);

	lin_approx_t us, uc;
	make_grid(us.x, 1, 2, 500000);
	uc.x = us.x;

	/*auto A = calcGlobalMatrix(us, c);
	auto b = calcB<EVector>(us, c, fs, fc);
	setFirstBoundaryConditions(A, b, us_true, uc_true, us);
	
	Eigen::JacobiSVD<EMatrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	EVector x = svd.solve(b);
	
	setAnswer(us, uc, x);*/

	auto A2 = calcGlobalMatrixProfile(us, c);
	auto b2 = to(calcB<Vector>(us, c, fs, fc));
	setFirstBoundaryConditions(A2, b2, us_true, uc_true, us);

	SLAU slau;
	slau.maxiter = 10;
	slau.eps = 1e-16;
	slau.is_log = true;
	slau.n = A2.n;
	slau.a = A2;
	slau.f = b2;
	slau.x.resize(slau.n);
	slau.t1.resize(slau.n);
	slau.t2.resize(slau.n);
	//slau.los2();
	slau.bsg_stab_lu();
	setAnswer(us, uc, to(slau.x));

	/*auto A1 = calcGlobalMatrixDiag(us, c);
	auto b1 = calcB<Vector>(us, c, fs, fc);
	setFirstBoundaryConditions(A1, b1, us_true, uc_true, us);*/
	/*Matrix denseA(A.cols(), A.rows());
	for (int i = 0; i < A.cols(); i++) {
		for (int j = 0; j < A.rows(); j++) {
			denseA(i, j) = A(i, j);
		}
	}
	MatrixDiagonal A2(denseA);*/

	/*Vector x1;
	SolverSLAE_Iterative solver;
	solver.epsilon = 1e-7;
	solver.isLog = false;
	solver.maxIterations = 5000;
	solver.start = Vector(us.size() * 2, 0);
	solver.w = 0.8;
	solver.seidel(A1, b1, x1);

	setAnswer(us, uc, x1);

	Vector b2;
	mul(A1, x1, b2);
	b2.negate();
	sum(b1, b2, b2);
	cout << "residual slae: " << calcNorm(b2) << endl;*/

	/*Matrix dense; A1.toDenseMatrix(dense);
	cout << A << endl << b << endl;
	dense.save(cout); cout << endl;
	b1.save(cout); cout << endl;
	cout << x << endl;
	x1.save(cout); cout << endl;*/

	/*Matrix dense; A2.toDense(dense);
	cout << A << endl;
	dense.save(cout); cout << endl;*/

	double residual = norm(us_true, us) + norm(uc_true, uc);

	/*cout << "answer s:    " << us.q << endl;
	cout << "should be s: " << calcTrulyApprox(us.x, us_true).q << endl;

	cout << "answer s:    " << uc.q << endl;
	cout << "should be s: " << calcTrulyApprox(uc.x, uc_true).q << endl;*/

	cout << "residual: " << residual << endl;

	system("pause");
}