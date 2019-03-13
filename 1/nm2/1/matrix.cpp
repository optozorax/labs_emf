#include <fstream>
#include <iomanip>
#include "matrix.h"

//-----------------------------------------------------------------------------
Matrix::Matrix(int n, int m, real fill) : m_matrix(m, std::vector<real>(n, fill)), m_n(n), m_m(m) {
}

//-----------------------------------------------------------------------------
void Matrix::loadFromFile(std::string fileName) {
	std::ifstream fin(fileName);

	m_matrix.clear();
	int n, m;
	fin >> n >> m;
	resize(n, m);
	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < width(); ++j) {
			fin >> operator()(i, j);
		}
	}

	fin.close();
}

//-----------------------------------------------------------------------------
void Matrix::saveToFile(std::string fileName) const {
	std::ofstream fout(fileName);

	fout.precision(std::numeric_limits<real>::digits10);
	int w = std::numeric_limits<real>::digits10 + 4;
	save(fout);

	fout.close();
}

//-----------------------------------------------------------------------------
void Matrix::load(std::istream& in) {
	m_matrix.clear(); m_m = 0; m_n = 0;
	int n, m;
	in >> n >> m;
	resize(n, m);
	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < width(); ++j) {
			in >> operator()(i, j);
		}
	}
}

//-----------------------------------------------------------------------------
void Matrix::save(std::ostream& out) const {
	out << m_n << "\t" << m_m << std::endl;
	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < width(); ++j)
			out << "\t" << std::setw(10) << operator()(i, j);
		out << std::endl;
	}
	out << std::endl;
}


//-----------------------------------------------------------------------------
void Matrix::getFromVector(int n, int m, const std::vector<real>& data) {
	resize(n, m);
	for (int i = 0; i < data.size(); ++i)
		operator()(i / n, i % n) = data[i];
}

//-----------------------------------------------------------------------------
void Matrix::resize(int n, int m, real fill) {
	if (m_n != n || m_m != m) {
		m_n = n;
		m_m = m;
		m_matrix.clear();
		m_matrix.resize(m_m, std::vector<real>(m_n, fill));
	}
}

//-----------------------------------------------------------------------------
void Matrix::negate(void) {
	for (auto& i : m_matrix) {
		for (auto& j : i) {
			j = -j;
		}
	}
}

//-----------------------------------------------------------------------------
bool Matrix::isSymmetric(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j <= i; ++j) {
			const real& a = operator()(i, j);
			const real& b = operator()(j, i);
			if (!isNear(a, b))
				return false;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isLowerTriangular(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < i; ++j) {
			if (fabs(operator()(j, i)) > 0.000001)
				return false;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isUpperTriangular(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < i; ++j) {
			if (operator()(i, j) != 0)
				return false;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isDiagonal(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		for (int j = 0; j < i; ++j) {
			if (operator()(j, i) != 0 && operator()(i, j) != 0)
				return false;
		}
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isDiagonalIdentity(void) const {
	if (height() != width())
		return false;

	for (int i = 0; i < height(); ++i) {
		if (operator()(i, i) != 1)
			return false;
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool Matrix::isDegenerate(void) const {
	// TODO
	return false;
}

//-----------------------------------------------------------------------------
real& Matrix::operator()(int i, int j) {
	return m_matrix[i][j];
}

//-----------------------------------------------------------------------------
const real& Matrix::operator()(int i, int j) const {
	return m_matrix[i][j];
}

//-----------------------------------------------------------------------------
int Matrix::width(void) const {
	return m_n;
}

//-----------------------------------------------------------------------------
int Matrix::height(void) const {
	return m_m;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void generateSparseSymmetricMatrix(int n, int min, int max, real percent, Matrix& result) {
	result.resize(n, n, 0);

	int count = percent * n * n;

	for (int k = 0; k < count; ++k) {
		int i = intRandom(0, n-1);
		int j = intRandom(0, i-1);
		result(i, j) = intRandom(min, max);
		result(j, i) = result(i, j);
	}

	for (int i = 0; i < n; i++)
		result(i, i) = intRandom(1, max - min);
}

//-----------------------------------------------------------------------------
void generateLMatrix(int n, int min, int max, real percent, Matrix& result) {
	result.resize(n, n, 0);

	int count = percent * n * n / 2;

	for (int k = 0; k < count; ++k) {
		int i = intRandom(0, n-1);
		int j = intRandom(0, i-1);
		result(i, j) = intRandom(min, max);
	}

	for (int i = 0; i < n; ++i) {
		result(i, i) = 1;
	}
}

//-----------------------------------------------------------------------------
void generateDiagonalMatrix(int n, int min, int max, Matrix& result) {
	result.resize(n, n, 0);

	for (int i = 0; i < n; ++i)
		result(i, i) = intRandom(min, max);
}

//-----------------------------------------------------------------------------
void generateVector(int n, int min, int max, Matrix& result) {
	result.resize(1, n, 0);

	for (int i = 0; i < n; ++i)
		result(i, 0) = intRandom(min, max);
}

//-----------------------------------------------------------------------------
void generateVector(int n, Matrix& result) {
	result.resize(1, n, 0);

	for (int i = 0; i < n; ++i)
		result(i, 0) = i+1;
}

//-----------------------------------------------------------------------------
void generateGilbertMatrix(int n, Matrix& result) {
	result.resize(n, n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			result(i, j) = double(1.0)/double((i+1)+(j+1)-1);
		}
	}
}

//-----------------------------------------------------------------------------
void generateTestMatrix(int n, int profileSize, Matrix& result) {
	result.resize(n, n);

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < profileSize; ++j) if (i-j-1 >= 0) {
			result(i, i-j-1) = -intRandom(0, 5);
			result(i-j-1, i) = result(i, i-j-1);
		}
	}

	for (int i = 0; i < n; ++i) {
		sumreal sum = 0;
		for (int j = 0; j < n; ++j) if (i != j) {
			sum += result(i, j);
		}
		result(i, i) = -sum;
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
bool mul(const Matrix& a, const Matrix& b, Matrix& result) {
	// result = rus_a * b
	if (a.width() != b.height())
		return false;

	result.resize(b.width(), a.height());

	for (int i = 0; i < b.width(); ++i) {
		for (int j = 0; j < a.height(); ++j) {
			real sum = 0;
			for (int k = 0; k < a.width(); ++k) {
				sum += a(j, k) * b(k, i);
			}
			result(j, i) = sum;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool sum(const Matrix& a, const Matrix& b, Matrix& result) {
	// result = rus_a + b
	if (a.width() != b.width() || a.height() != b.height())
		return false;

	result.resize(a.width(), a.height());

	for (int i = 0; i < a.width(); ++i) {
		for (int j = 0; j < a.height(); ++j) {
			result(j, i) = a(j, i) + b(j, i);
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
bool transpose(Matrix& a) {
	// rus_a = rus_a^T
	Matrix a_t(a.height(), a.width());

	for (int i = 0; i < a.height(); ++i) {
		for (int j = 0; j < a.width(); ++j) {
			a_t(j, i) = a(i, j);
		}
	}

	a = a_t;

	return true;
}

//-----------------------------------------------------------------------------
sumreal sumAllElementsAbs(const Matrix& a) {
	sumreal sum = 0;
	for (int i = 0; i < a.height(); ++i) {
		for (int j = 0; j < a.width(); ++j) {
			sum += fabs(a(i, j));
		}
	}

	return sum;
}

//-----------------------------------------------------------------------------
bool calcLUsq(const Matrix& a, Matrix& l, Matrix& u) {
	if (a.width() != a.height())
		return false;

	l.resize(a.width(), a.height(), 0);
	u.resize(a.width(), a.height(), 0);	

	for (int i = 0; i < a.height(); ++i) {
		// Считаем элементы матрицы L
		for (int j = 0; j < i; ++j) {
			real sum = 0;
			for (int k = 0; k < j; ++k)
				sum += l(i, k) * u(k, j);

			l(i, j) = (a(i, j) - sum) / l(j, j);
		}

		// Считаем элементы матрицы U
		for (int j = 0; j < i; ++j) {
			real sum = 0;
			for (int k = 0; k < j; ++k)
				sum += l(j, k) * u(k, i);

			u(j, i) = (a(j, i) - sum) / u(j, j);
		}

		// Считаем диагональный элемент
		real sum = 0;
		for (int k = 0; k < i; ++k)
			sum += l(i, k) * u(k, i);

		l(i, i) = sqrt(a(i, i) - sum);
		u(i, i) = l(i, i);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool calcLUsq_partial(const Matrix& a, Matrix& l, Matrix& u) {
	if (a.width() != a.height())
		return false;

	l.resize(a.width(), a.height(), 0);
	u.resize(a.width(), a.height(), 0);	

	for (int i = 0; i < a.height(); ++i) {
		// Считаем элементы матрицы L
		for (int j = 0; j < i; ++j) if (a(i, j) != 0) {
			real sum = 0;
			for (int k = 0; k < j; ++k)
				sum += l(i, k) * u(k, j);

			l(i, j) = (a(i, j) - sum) / l(j, j);
		}

		// Считаем элементы матрицы U
		for (int j = 0; j < i; ++j) if (a(j, i) != 0) {
			real sum = 0;
			for (int k = 0; k < j; ++k)
				sum += l(j, k) * u(k, i);

			u(j, i) = (a(j, i) - sum) / u(j, j);
		}

		// Считаем диагональный элемент
		real sum = 0;
		for (int k = 0; k < i; ++k)
			sum += l(i, k) * u(k, i);

		l(i, i) = sqrt(a(i, i) - sum);
		u(i, i) = l(i, i);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
bool calcLDL(const Matrix& a, Matrix& l, Matrix& d) {
	// l * d * l^T = rus_a
	if (!a.isSymmetric())
		return false;

	l.resize(a.width(), a.height(), 0);
	d.resize(a.width(), a.height(), 0);	

	for (int i = 0; i < a.height(); ++i) {
		// Считаем элементы матрицы L
		for (int j = 0; j < i; ++j) {
			real sum = 0;
			for (int k = 0; k < j; ++k)
				sum += d(k, k) * l(j, k) * l(i, k);

			if (fabs(d(j, j)) < 0.0001)
				l(i, j) = 0;
			else
				l(i, j) = (a(i, j) - sum) / d(j, j);
		}

		// Считаем диагональный элемент
		{
			real sum = 0;
			for (int j = 0; j < i; ++j)
				sum += d(j, j) * l(i, j) * l(i, j);
			d(i, i) = a(i, i) - sum;
		}
	}

	for (int i = 0; i < l.height(); i++)
		l(i, i) = 1;
	
	return true;
}

//-----------------------------------------------------------------------------
bool calcGaussianReverseOrder(const Matrix& l, const Matrix& y, Matrix& x) {
	// l * x = y, l - нижнетреугольная матрица
	if (!l.isLowerTriangular() || !l.isDiagonalIdentity())	
		return false;

	x.resize(1, y.height());

	for (int i = x.height() - 1; i >= 0; --i) {
		real sum = 0;
		for (int j = i; j < x.height(); ++j)
			sum += l(j, i) * x(j, 0);
		x(i, 0) = y(i, 0) - sum;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool calcGaussianFrontOrder(const Matrix& l, const Matrix& y, Matrix& x) {
	// l * x = y, l - верхнетреугольная матрица
	if (!l.isLowerTriangular() || !l.isDiagonalIdentity())	
		return false;

	x.resize(1, y.height());

	for (int i = 0; i < x.height(); ++i) {
		real sum = 0;
		for (int j = 0; j < i; ++j)
			sum += l(i, j) * x(j, 0);
		x(i, 0) = y(i, 0) - sum;
	}

	return true;
}

//-----------------------------------------------------------------------------
bool calcGaussianCentralOrder(const Matrix& d, const Matrix& y, Matrix& x) {
	// d * x = y, d - диагональная матрица
	if (!d.isDiagonal())	
		return false;

	x.resize(1, y.height());

	for (int i = 0; i < x.height(); ++i)
		x(i, 0) = y(i, 0) / d(i, i);

	return true;
}

//-----------------------------------------------------------------------------
bool solveSLAE_by_LDL(const Matrix& a, const Matrix& y, Matrix& x) {
	// rus_a * x = y, rus_a - симметричная матрица
	if (!(a.width() == a.height() && a.width() == y.height() && !a.isDegenerate()))
		return false;

	Matrix l, d, z, w;

	if (!calcLDL(a, l, d))
		return false;

	if (!calcGaussianFrontOrder(l, y, z))
		return false;

	if (!calcGaussianCentralOrder(d, z, w))
		return false;

	if (!calcGaussianReverseOrder(l, w, x))
		return false;

	return true;
}

//-----------------------------------------------------------------------------
bool solveSLAE_byGaussMethod(const Matrix& a1, const Matrix& y1, Matrix& x1) {
	if (!(a1.width() == a1.height() && y1.height() == a1.width() && !a1.isDegenerate()))
		return false;

	Matrix a(a1);
	Matrix y(y1);

	for (int i = 0; i < a.height(); ++i) {
		// Находим максимальный элемент
		int maxI = i;
		for (int j = i+1; j < a.height(); ++j)
			if (fabs(a(j, i)) > fabs(a(maxI, i)))
				maxI = j;

		// Переставляем эту строчку с текущей
		for (int j = i; j < a.width(); ++j)
			std::swap(a(i, j), a(maxI, j));
		std::swap(y(i, 0), y(maxI, 0));

		// Перебираем все строчки ниже и отнимаем текущую строчку от них
		for (int j = i+1; j < a.height(); ++j) {
			real m = a(j, i) / a(i, i);
			for (int k = i; k < a.width(); ++k)
				a(j, k) -= m * a(i, k);
			y(j, 0) -= m * y(i, 0);
		}

		// Делим текущую строку на ее ведущий элемент, чтобы на диагонали были единицы
		double m = a(i, i);
		for (int j = i; j < a.width(); ++j)
			a(i, j) /= m;
		y(i, 0) /= m;
	}

	// Считаем обратный ход Гаусса
	transpose(a);
	calcGaussianReverseOrder(a, y, x1);

	return true;
}