#pragma once

/** Файл для работы с матрицей в разряженном формате и решении */

#include <vector>
#include <map>
#include <iostream>
#include "lib.h"

//-----------------------------------------------------------------------------
/** Класс квадратной разряженной матрицы в строчно-столбцовом формате с симметричным профилем. Примечание: ненулевым считается элемент, который имеется в профиле, неважно что в массивах l, u он может иметь значение 0. */
class matrix_sparse_t
{
public:
	int n;            /// Размерность матрицы
	vector<double> d; /// Диагональные элементы матрицы
	vector<double> l; /// Элементы матрицы из нижнего треугольника
	vector<double> u; /// Элементы матрицы из верхнего треугольника
	vector<int> i;    /// Массив начала строк в формате (ia в методичке)
	vector<int> j;    /// Массив столбцов каждого элемента (ja в методичке)

	matrix_sparse_t(int n);

	/** Преобразование разреженной матрицы к плотному формату. */
	void to_dense(matrix_t& m) const;

	void clear_line(int line);

	int line_elem_start(int line) const; /// Получить позицию в массивах l, u элемента, с которого начинается строка line
	int line_elem_row(int line, int elem) const; /// Получить столбец ненулевого элемента в строке line под номером elem
	int line_elem_count(int line) const; /// Получить количество ненулевых элементов в строке

	/** Раскладывает текущую матрицу в неполное LU разложение и хранит результат в матрице lu. Неполное разложение - это когда были применены формулы для получения LU матрицы, но только к существующим ненулевым элементам, без перестройки формата. Иными словами называется "неполная факторизация". */
	void decompose_lu_partial(matrix_sparse_t& lu) const;

	/* Методы для умножения разряженной матрицы на вектор. */
	void mul(vector_t& x_y) const;         // x_y = a * x_y
	void mul_t(vector_t& x_y) const;       // x_y = a^t * x_y

	/* Представляет, что текущая матрица хранит LU разложение, и соответственно можно каждую матрицу этого разложения умножить на соответствующие вектора. */
	void mul_l_inv_t(vector_t& x_y) const; // x_y = l^-t * x_y
	void mul_u_inv_t(vector_t& x_y) const; // x_y = u^-t * x_y
	void mul_l_inv(vector_t& x_y) const;   // x_y = l^-1 * x_y
	void mul_u_inv(vector_t& x_y) const;   // x_y = u^-1 * x_y
	void mul_u(vector_t& x_y) const;       // x_y = u * x_y	
};

ostream& operator<<(ostream& out, const matrix_sparse_t& m);

//-----------------------------------------------------------------------------
/** Квадратная матрица с произвольным доступом к любому элементу. Предполагается, что матрица будет разряженная. Далее можно перегенерировать её в разряженную матрицу. */
class matrix_sparse_ra_t
{
public:
	matrix_sparse_ra_t(int n);

	/** Установить значение в позиции (i, j) */
	double& operator()(int i, int j);

	/** Получить значение в позиции (i, j). Если туда ещё не устанавливалось значение, вызывается исключение. */
	const double& operator()(int i, int j) const;

	/** Преобразует текущую матрицу к разреженной матрице. */
	matrix_sparse_t to_sparse(void) const;
private:
	int n;
	vector<double> dm;
	vector<map<int, double>> lm, um;
};

//-----------------------------------------------------------------------------
/* Функции для умножения векторов. */
void mul(const vector_t& d, vector_t& x_y);     // x_y = d * x_y
void mul_inv(const vector_t& d, vector_t& x_y); // x_y = x_y / d

//-----------------------------------------------------------------------------
/** Решает СЛАУ с матрицей в разряженном формате при помощи Локально-Оптимальной Схемы (ЛОС) с предобуславливанием на основе неполной LU факторизации. */
vector_t solve_by_los_lu(
	const matrix_sparse_t& a,
	const vector_t& b,
	int maxiter,
	double eps,
	bool is_log = false
);