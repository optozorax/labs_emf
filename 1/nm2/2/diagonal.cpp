#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "diagonal.h"

//-----------------------------------------------------------------------------
Diagonal::Diagonal(int n) : n(n) {
}

//-----------------------------------------------------------------------------
int Diagonal::calcDiagonalsCount(void) {
	return 2 * n - 1;
}

//-----------------------------------------------------------------------------
int Diagonal::calcMinDiagonal(void) {
	return -(n-1);
}

//-----------------------------------------------------------------------------
int Diagonal::calcMaxDiagonal(void) {
	return n - 1;
}

//-----------------------------------------------------------------------------
int Diagonal::calcDiagonalSize(int d) {
	return n - std::abs(d);
}

//-----------------------------------------------------------------------------
bool Diagonal::isLineIntersectDiagonal(int line, int d) {
	if (d <= 0)
		return (line+d >= 0);

	if (d > 0)
		return (line < calcDiagonalSize(d));
}

//-----------------------------------------------------------------------------
bool Diagonal::isRowIntersectDiagonal(int row, int d) {
	return isLineIntersectDiagonal(row, -d);
}

//-----------------------------------------------------------------------------
int Diagonal::calcLine_byDP(int d, int pos) {
	if (d <= 0)
		return -d + pos;

	if (d > 0)
		return pos;
}

//-----------------------------------------------------------------------------
int Diagonal::calcRow_byDP(int d, int pos) {
	if (d <= 0)
		return pos;

	if (d > 0)
		return pos + d;
}

//-----------------------------------------------------------------------------
int Diagonal::calcDiag_byLR(int line, int row) {
	return row - line;
}

//-----------------------------------------------------------------------------
int Diagonal::calcPos_byLR(int line, int row) {
	return calcPos_byDL(calcDiag_byLR(line, row), line);
}

//-----------------------------------------------------------------------------
int Diagonal::calcPos_byDL(int d, int line) {
	if (d <= 0)
		return line+d;

	if (d > 0)
		return line;
}

//-----------------------------------------------------------------------------
int Diagonal::calcPos_byDR(int d, int row) {
	return calcPos_byDL(d, calcLine_byDR(d, row));
}

//-----------------------------------------------------------------------------
int Diagonal::calcRow_byDL(int d, int line) {
	return line+d;
}

//-----------------------------------------------------------------------------
int Diagonal::calcLine_byDR(int d, int row) {
	return calcRow_byDL(-d, row);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
MatrixDiagonal::MatrixDiagonal() : n(0), dc(n) {
}

//-----------------------------------------------------------------------------
MatrixDiagonal::MatrixDiagonal(int n, std::vector<int> format) : dc(n) {
	resize(n, format);
}

//-----------------------------------------------------------------------------
MatrixDiagonal::MatrixDiagonal(const Matrix& a) : dc(n) {
	if (a.width() != a.height())
		throw std::exception();

	n = a.width();
	dc.n = n;

	// Определяем формат
	std::vector<int> format;
	format.clear();
	format.push_back(0);
	for (int i = dc.calcMinDiagonal(); i <= dc.calcMaxDiagonal(); ++i) 
		if (i != 0) {
			auto mit = matrix_diagonal_iterator(n, i, false);
			auto mite = matrix_diagonal_iterator(n, i, true);
			for (; mit != mite; ++mit) {
				if (a(mit.i, mit.j) != 0) {
					format.push_back(i);
					break;
				}
			}
		}

	// Создаем формат
	resize(n, format);

	// Обходим массив и записываем элементы
	for (int i = 0; i < getDiagonalsCount(); ++i) {
		auto mit = posBegin(i);
		for (auto it = begin(i); it != end(i); ++it, ++mit)
			*it = a(mit.i, mit.j);
	}
}

//-----------------------------------------------------------------------------
void MatrixDiagonal::toDenseMatrix(Matrix& dense) const {
	dense.resize(n, n, 0);

	// Обходим массив и записываем элементы
	for (int i = 0; i < getDiagonalsCount(); ++i) {
		auto mit = posBegin(i);
		for (auto it = begin(i); it != end(i); ++it, ++mit)
			dense(mit.i, mit.j) = *it;
	}
}

//-----------------------------------------------------------------------------
void MatrixDiagonal::resize(int n1, std::vector<int> format) {
	if (format[0] != 0)
		throw std::exception();

	dc.n = n1;
	n = n1;
	fi = format;

	di.clear();
	for (const auto& i : format)
		di.push_back(std::vector<real>(dc.calcDiagonalSize(i), 0));
}

//-----------------------------------------------------------------------------
void MatrixDiagonal::save(std::ostream& out) const {
	out << n << " " << fi.size() << std::endl;
	for (const auto& i : fi)
		out << i << " ";
	out << std::endl;

	for (int i = 0; i < getDiagonalsCount(); ++i) {
		for (auto j = begin(i); j != end(i); ++j) {
			out << (*j) << " ";
		}
		out << std::endl;
	}
	out << std::endl;
}

//-----------------------------------------------------------------------------
void MatrixDiagonal::load(std::istream& in) {
	int n, m;
	in >> n >> m;
	std::vector<int> format(m, 0);
	for (int i = 0; i < m; ++i)
		in >> format[i];
	resize(n, format);
	for (int i = 0; i < getDiagonalsCount(); ++i) {
		for (auto j = begin(i); j != end(i); ++j) {
			in >> (*j);
		}
	}
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::dimension(void) const {
	return n;
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::getDiagonalsCount(void) const {
	return di.size();
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::getDiagonalSize(int diagNo) const {
	return di[diagNo].size();
}

//-----------------------------------------------------------------------------
int MatrixDiagonal::getDiagonalPos(int diagNo) const {
	return fi[diagNo];
}

//-----------------------------------------------------------------------------
std::vector<int> MatrixDiagonal::getFormat(void) const {
	return fi;
}

//-----------------------------------------------------------------------------
matrix_diagonal_iterator MatrixDiagonal::posBegin(int diagNo) const {
	return matrix_diagonal_iterator(n, fi[diagNo], false);
} 

//-----------------------------------------------------------------------------
matrix_diagonal_iterator MatrixDiagonal::posEnd(int diagNo) const {
	return matrix_diagonal_iterator(n, fi[diagNo], true);
}

//-----------------------------------------------------------------------------
MatrixDiagonal::iterator MatrixDiagonal::begin(int diagNo) {
	return di[diagNo].begin();
}

//-----------------------------------------------------------------------------
MatrixDiagonal::const_iterator MatrixDiagonal::begin(int diagNo) const {
	return di[diagNo].begin();
}

//-----------------------------------------------------------------------------
MatrixDiagonal::iterator MatrixDiagonal::end(int diagNo) {
	return di[diagNo].end();
}

//-----------------------------------------------------------------------------
MatrixDiagonal::const_iterator MatrixDiagonal::end(int diagNo) const {
	return di[diagNo].end();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
matrix_diagonal_iterator::matrix_diagonal_iterator(int n, int d, bool isEnd) {
	Diagonal dc(n);
	if (isEnd) {
		i = dc.calcLine_byDP(d, dc.calcDiagonalSize(d));
		j = dc.calcRow_byDP(d, dc.calcDiagonalSize(d));
	} else {
		i = dc.calcLine_byDP(d, 0);
		j = dc.calcRow_byDP(d, 0);
	}
}

//-----------------------------------------------------------------------------
matrix_diagonal_iterator& matrix_diagonal_iterator::operator++() {
	i++;
	j++;
	return *this;
}

//-----------------------------------------------------------------------------
matrix_diagonal_iterator matrix_diagonal_iterator::operator++(int) {
	i++;
	j++;
	return *this;
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_iterator::operator==(const matrix_diagonal_iterator& b) const {
	return b.i == i && b.j == j;
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_iterator::operator!=(const matrix_diagonal_iterator& b) const {
	return b.i != i || b.j != j;
}

//-----------------------------------------------------------------------------
matrix_diagonal_iterator& matrix_diagonal_iterator::operator+=(const ptrdiff_t& movement) {
	i += movement;
	j += movement;
	return *this;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
matrix_diagonal_line_iterator::matrix_diagonal_line_iterator(int n, std::vector<int> format, bool isOnlyLowTriangle) : dc(n), m_isEnd(false), m_isLineEnd(false) {
	// Создаем обратное преобразование из формата диагонали в ее номер в формате
	for (int i = 0; i < format.size(); ++i)
		if ((isOnlyLowTriangle && format[i] < 0) || !isOnlyLowTriangle)
			m_map[format[i]] = i;

	// Создаем сортированный формат, чтобы по нему двигаться
	if (isOnlyLowTriangle) {
		for (int i = 0; i < format.size(); ++i)
			if (format[i] < 0)
				m_sorted_format.push_back(format[i]);
	} else 
		m_sorted_format = format;
	std::sort(m_sorted_format.begin(), m_sorted_format.end());

	line = 0;
	pos = 0;
	start = m_sorted_format.size() - 1;
	end = start;

	// Находим, с какой диагонали начинается текущая строка
	for (int i = 0; i < m_sorted_format.size(); ++i) {
		if (dc.isLineIntersectDiagonal(line, m_sorted_format[i])) {
			start = i;
			break;
		}
	}

	// Находим на какой диагонали кончается текущая строка
	for (int i = 0; i < m_sorted_format.size(); ++i) {
		int j = m_sorted_format.size() - i - 1;
		if (dc.isLineIntersectDiagonal(line, m_sorted_format[j])) {
			end = j;
			break;
		}
	}

	calcPos();
}

//-----------------------------------------------------------------------------
matrix_diagonal_line_iterator& matrix_diagonal_line_iterator::operator++() {
	if (!m_isEnd) {
		if (m_isLineEnd) {
			// Сдвигаемся на одну строку
			line++;

			// Определяем какие диагонали пересекают эту строку
			if (start != 0)
				if (dc.isLineIntersectDiagonal(line, m_sorted_format[start-1]))
					start = start-1;

			if (end != 0)
				if (!dc.isLineIntersectDiagonal(line, m_sorted_format[end]))
					if (start != end)
						end = end-1;

			m_isLineEnd = false;
			if (line == dc.n)
				m_isEnd = true;

			pos = 0;
			calcPos();
		} else {
			// Сдвигаемся на один столбец
			pos++;
			calcPos();
		}
	}

	return *this;
}

//-----------------------------------------------------------------------------
matrix_diagonal_line_iterator matrix_diagonal_line_iterator::operator++(int) {
	return operator++();
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_line_iterator::isLineEnd(void) const {
	return m_isLineEnd;
}

//-----------------------------------------------------------------------------
bool matrix_diagonal_line_iterator::isEnd(void) const {
	return m_isEnd;
}

//-----------------------------------------------------------------------------
void matrix_diagonal_line_iterator::calcPos(void) {
	// Вычисляет все текущие положения согласно переменным start, pos и формату
	if (!dc.isLineIntersectDiagonal(line, m_sorted_format[end]) || (start + pos > end)) {
		m_isLineEnd = true;
		i = line;
		j = 0;
		d = 0;
		di = 0;
		dn = 0;
	} else {
		i = line;
		d = m_sorted_format[start + pos];
		dn = m_map[d];
		di = dc.calcPos_byDL(d, i);
		j = dc.calcRow_byDL(d, i);
	}
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
std::vector<int> makeSevenDiagonalFormat(int n, int m, int k) {
	std::vector<int> result;

	if (1+m+k >= n)
		throw std::exception();

	result.push_back(0);

	result.push_back(1);
	result.push_back(1+m);
	result.push_back(1+m+k);

	result.push_back(-1);
	result.push_back(-1-m);
	result.push_back(-1-m-k);

	return result;
}

//-----------------------------------------------------------------------------
std::vector<int> generateRandomFormat(int n, int diagonalsCount) {
	Diagonal d(n);

	std::vector<int> result;
	result.push_back(0);

	// Создаем массив всех возможных диагоналей
	std::vector<int> diagonals;
	for (int i = d.calcMinDiagonal(); i <= d.calcMaxDiagonal(); ++i)
		if (i != 0)
			diagonals.push_back(i);

	diagonalsCount = std::min<int>(diagonals.size(), diagonalsCount);

	// Заполняем результат случайными диагоналями из этого массива
	for (int i = 0; i < diagonalsCount; ++i) {
		int pos = intRandom(0, diagonals.size());
		result.push_back(diagonals[pos]);
		diagonals.erase(diagonals.begin() + pos);
	}

	return result;
}

//-----------------------------------------------------------------------------
void generateDiagonalMatrix(int n, int min, int max, std::vector<int> format, MatrixDiagonal& result) {
	result.resize(n, format);
	for (int i = 0; i < result.getDiagonalsCount(); ++i) {
		auto mit = result.posBegin(i);
		for (auto it = result.begin(i); it != result.end(i); ++it, ++mit)
			(*it) = intRandom(min, max);
	}
}

//-----------------------------------------------------------------------------
void generateDiagonallyDominantMatrix(int n, std::vector<int> format, bool isNegative,  MatrixDiagonal& result) {
	result.resize(n, format);

	for (int i = 0; i < result.getDiagonalsCount(); ++i) {
		auto mit = result.posBegin(i);
		for (auto it = result.begin(i); it != result.end(i); ++it, ++mit) {
			if (isNegative)
				*it = -intRandom(0, 5);
			else
				*it = intRandom(0, 5);
		}
	}

	matrix_diagonal_line_iterator mit(n, format, false);
	for (; !mit.isEnd(); ++mit) {
		sumreal& sum = result.begin(0)[mit.i];
		sum = 0;
		for (; !mit.isLineEnd(); ++mit)
			if (mit.i != mit.j)
				sum += result.begin(mit.dn)[mit.di];
		sum = std::fabs(sum);
	}

	result.begin(0)[0] += 1;
}

//-----------------------------------------------------------------------------
bool mul(const MatrixDiagonal& a, const Vector& x, Vector& y) {
	if (x.size() != a.dimension())
		return false;

	y.resize(x.size());

	// Зануление результата
	y.zero();
	for (int i = 0; i < a.getDiagonalsCount(); ++i) {
		auto mit = a.posBegin(i);
		for (auto it = a.begin(i); it != a.end(i); ++it, ++mit)
			y(mit.i) += (*it) * x(mit.j);
	}

	return true;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
SolverSLAE_Iterative::SolverSLAE_Iterative() : 
	w(1), 
	isLog(false), 
	log(std::cout), 
	start(), 
	epsilon(0.00001), 
	maxIterations(100) {
}

//-----------------------------------------------------------------------------
void SolverSLAE_Iterative::save(std::ostream& out) const {
	out << w << std::endl;
	out << isLog << std::endl;
	start.save(out);
	out << std::scientific;
	out << epsilon << std::endl;
	out << std::defaultfloat;
	out << maxIterations << std::endl;
}

//-----------------------------------------------------------------------------
void SolverSLAE_Iterative::load(std::istream& in) {
	in >> w >> isLog;
	start.load(in);
	in >> epsilon >> maxIterations;
}

//-----------------------------------------------------------------------------
IterationsResult SolverSLAE_Iterative::jacobi(const MatrixDiagonal& a, const Vector& y, Vector& x) const {
	return iteration_process(a, y, x, &SolverSLAE_Iterative::iteration_jacobi);
}

//-----------------------------------------------------------------------------
IterationsResult SolverSLAE_Iterative::seidel(const MatrixDiagonal& a, const Vector& y, Vector& x) const {
	return iteration_process(a, y, x, &SolverSLAE_Iterative::iteration_seidel);
}

//-----------------------------------------------------------------------------
IterationsResult SolverSLAE_Iterative::iteration_process(const MatrixDiagonal& a, const Vector& y, Vector& x, step_function step) const {
	if (a.dimension() != y.size() || start.size() != y.size())
		throw std::exception();

	// Считаем норму матрицы: ее максимальный элемент по модулю
	real yNorm = calcNorm(y);
	x1.resize(y.size());
	x = start;

	// Цикл по итерациям
	int i = 0;
	real relativeResidual = epsilon + 1;
	for (; i < maxIterations && relativeResidual > epsilon; ++i) {
		// Итерационный шаг
		step(this, a, y, x);

		// Считаем невязку
		mul(a, x, x1);
		x1.negate();
		sum(x1, y, x1);
		relativeResidual = fabs(calcNorm(x1)) / yNorm;

		// Выводим данные
		if (isLog)
			log << i << "\t" << std::scientific << std::setprecision(3) << relativeResidual << std::endl;
	}

	return {i, relativeResidual};
}

//-----------------------------------------------------------------------------
void SolverSLAE_Iterative::iteration_jacobi(const MatrixDiagonal& a, const Vector& y, Vector& x) const {
	// Умножаем матрицу на решение
	mul(a, x, x1);

	// x^(k+1) = x^k + w/a(i, i) * x^(k+1)
	auto it = a.begin(0);
	for (int i = 0; i < x1.size(); ++i, ++it)
		x(i) += w / (*it) * (y(i)-x1(i));
}

//-----------------------------------------------------------------------------
void SolverSLAE_Iterative::iteration_seidel(const MatrixDiagonal& a, const Vector& y, Vector& x) const {
	// Умножем верхний треугольник на решение
	x1.zero();
	for (int i = 0; i < a.getDiagonalsCount(); ++i)
		if (a.getDiagonalPos(i) >= 0) {
			auto mit = a.posBegin(i);
			for (auto it = a.begin(i); it != a.end(i); ++it, ++mit)
				x1(mit.i) += (*it) * x(mit.j);
		}

	// Проходим по нижнему треугольнику и считаем все параметры
	matrix_diagonal_line_iterator mit(a.dimension(), a.getFormat(), true);
	for (; !mit.isEnd(); ++mit) {
		for (; !mit.isLineEnd(); ++mit)
			x1(mit.i) += a.begin(mit.dn)[mit.di] * x(mit.j);
		x(mit.i) = x(mit.i) + w/a.begin(0)[mit.i] * (y(mit.i) - x1(mit.i));
	}
}