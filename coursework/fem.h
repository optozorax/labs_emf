#pragma once

/* Заголовок функций для реализации Метода Конечных Элементов (МКЭ) в 2D пространстве. Уравнение гиперболическое (с второй производной по времени). Схема для аппроксимации по времени: Кранка-Николсона. Базисные элементы: билинейные. Форма сетки: прямоугольники. */

/* Сделано только для того, чтобы сдать курсовую, поэтому никакой расширяемости и мега-красивого кода здесь ожидать не стоит. Разве что можете использовать в своих курсовых, я старался писать как можно более красивей с учетом того, что я начал делать курсовую в воскресенье, а сдавать мне её в среду. */

#include "lib.h"
#include "sparse.h"

//-----------------------------------------------------------------------------
/** Узел конечного элемента. Другими словами, вес, домноженный на базовую функцию. Из сумм этих элементов образуется конечный элемент. */
struct basic_elem_t
{
	int i; /// Номер узла

	double x, y; /// Координаты узла

	basic_elem_t *up, *down, *left, *right; /// Указатели на соседей узла

	/** Проверяет, является ли элемент граничным. Он таким явлется, если у него нет хотя бы одного соседа. */
	bool is_boundary(void) const;
};

/** Прямоугольный конечный элемент на основе билинейных базисных функций. Образуется из четырех узлов. */
struct elem_t
{
	int i; /// Номер конечного элемента
	basic_elem_t* e[4]; /// Указатели на все 4 элемента конечного узла, нумерация такая:
	/** 
		 Y
		 ^  3 +-----+ 4
		 |    |     |
		 |  1 +-----+ 2
		-+--------------> X
		 |
	 */

	double get_hx(void) const; /// Ширина конечного элемента
	double get_hy(void) const; /// Высота конечного элемента

	/** Рассчитать значение внутри конечного элемента. q - вектор рассчитыванных весов. */
	double value(double x, double y, const vector_t& q) const;
};

/** Все константы решаемого уравнения. */
struct constants_t
{
	double lambda; /// Коэффициент внутри div
	double gamma;  /// Коэффициент при u
	double sigma;  /// Коэффициент при du/dt
	double chi;    /// Коэффициент при d^2u/dt^2
};

//-----------------------------------------------------------------------------
/** t in [-1, 1]. x in [0, 1] При t=-1 возвращаемое значение полностью смещается к 0, при t=1 возвращаемое значение полностью смещается к 1, при t=0 возвращаемое значение равно x. Между этими значениями используются формулы, чтобы решение сгущалось постепенно к одному из концов. Функция используется для генерации неравномернной сетки. */
double non_linear_offset(double x, double t);

/** Генерирует неравномерную сетку по заданным параметрам. n - число внутренних узлов. То есть если n будет равно 0, то узел под номером 0 будет a, а под номером 1 будет b. */
class grid_generator_t
{
public:
	grid_generator_t(double a, double b, int n, double t = 0);
	double operator()(int i) const;
	int size() const;
private:
	double a, len, t, n1;
};

/** Класс двумерной неравномерной сетки по пространству в виде прямоугольника. */
class grid_t
{
public:
	vector<elem_t> es; /// Массив конечных элементов сетки
	vector<basic_elem_t> bes; /// Массив узлов сетки
	int n; /// Число узлов

	/** Рассчитать неравномерную сетку. */
	void calc(const grid_generator_t& gx, const grid_generator_t& gy);
};

/** Рассчитать веса идеальной аппроксимации функции u при помощи узлов bes. */
vector_t calc_true_approx(const function_2d_t& u, const vector<basic_elem_t>& bes);

//-----------------------------------------------------------------------------
/* Расчет локальных матриц для конечного элемента. */
matrix_t calc_local_matrix_g(const elem_t& e, const constants_t& cs);
matrix_t calc_local_matrix_c(const elem_t& e);
matrix_t calc_local_matrix_m(const elem_t& e, const constants_t& cs);
vector_t calc_local_vector_b(const elem_t& e, const function_2d_t& f);

//-----------------------------------------------------------------------------
/** Рассчитать глобальный вектор из локальных векторов для всех конечных элементов. */
vector_t calc_global_vector(
	const vector<elem_t>& es,
	const function<vector_t(const elem_t&)> calc_local_vector,
	int n
);

/** Рассчитать глобальную матрицу из функции построения локальных матриц. */
matrix_sparse_t calc_global_matrix(
	const vector<elem_t>& es,
	const function<matrix_t(const elem_t&)> calc_local_matrix,
	int n
);

//-----------------------------------------------------------------------------
/* Численный расчет определенных интегралов. */
double calc_integral_gauss3(
	double a, double b, int n, // n - количество внутренных узлов
	const function_1d_t& f
);
double calc_integral_gauss3(
	double ax, double bx, int nx,
	double ay, double by, int ny,
	const function_2d_t& f
);

//-----------------------------------------------------------------------------
/* Численный расчет производной. */
function_1d_t calc_first_derivative(const function_1d_t& f);
function_1d_t calc_second_derivative(const function_1d_t& f);

//-----------------------------------------------------------------------------
/** Для гиперболического дифференциального уравнения и функции u считает каким должно быть f, чтобы решением этого диф. уравнения была функци u. Делает это численно. */
function_3d_t calc_right_part(const function_3d_t& u, const constants_t& cs);

//-----------------------------------------------------------------------------
/** Использует схему Кранка-Николсона для получения СЛАУ. Предполагается, что разряженные матрицы имеют одинаковый формат. */
void calc_crank_nicolson_method(
	const matrix_sparse_t& c,
	const matrix_sparse_t& m,
	const matrix_sparse_t& g,
	const vector_t& b0,  // b current (b_j)
	const vector_t& bl,  // b last (b_{j-1})
	const vector_t& bll, // b last last (b_{j-2})
	const vector_t& ql,  // q last (q_{j-1})
	const vector_t& qll, // q last last (q_{j-2})
	const constants_t& cs,
	const grid_generator_t& time_grid,
	int time_i,
	matrix_sparse_t& a,
	vector_t& b
);

//-----------------------------------------------------------------------------
/** Функция, которая устанавливает краевые условия для задачи в СЛАУ. Сделана для того, чтобы не посылать в функцию решения МКЭ истинную функцию, а чтобы посылать красивую оболочку, которую потенциально можно использовать в реальных задачах. */
typedef function<void(matrix_sparse_t&, vector_t&, const vector<basic_elem_t>&, double)> boundary_setter_t;

/** Записывает первые краевые условия в матрицу a и вектор b. Для этой записи ему необходимо получить истинную функцию. */
void write_first_boundary_conditions(
	matrix_sparse_t& a,
	vector_t& b,
	const vector<basic_elem_t>& bes,
	double t,
	const function_3d_t& u
);

//-----------------------------------------------------------------------------
/** Решает при помощи МКЭ дифференциальное уравнение с функцией правой части f, заданными константами, прямоугольной сеткой grid и функцией выставления краевых условий. Использует схему Кранка-Николсона для аппроксимации по времени, и ЛОС в разряженной строчно-столбцовой матрице для решения СЛАУ. */
vector<vector_t> solve_differential_equation(
	const function_3d_t& f,
	const boundary_setter_t& set_boundary_conditions,
	const vector_t& q0,
	const vector_t& q1,
	const constants_t& cs,
	const grid_t& grid,
	const grid_generator_t& time_grid
);