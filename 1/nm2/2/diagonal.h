#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <functional>
#include "../1/common.h"
#include "../1/vector.h"
#include "../1/matrix.h"

class MatrixDiagonal;
class matrix_diagonal_iterator;
class SolverSLAE_Iterative;

//-----------------------------------------------------------------------------
/** Класс для вычислений различных параметров с диагональными матрицами. */
class Diagonal
{
public:
	int n;

	//-------------------------------------------------------------------------
	Diagonal(int n);

	int calcDiagonalsCount(void);
	int calcMinDiagonal(void);
	int calcMaxDiagonal(void);
	int calcDiagonalSize(int d);

	bool isLineIntersectDiagonal(int line, int d);
	bool isRowIntersectDiagonal(int row, int d);

	//-------------------------------------------------------------------------
	/*
		R - Row - столбец
		L - Line - строка
		P - Pos - номер элемента в диагонали
		D - Diag - формат диагонали
	*/

	int calcLine_byDP(int d, int pos);
	int calcRow_byDP(int d, int pos);

	int calcDiag_byLR(int line, int row);
	int calcPos_byLR(int line, int row);

	int calcPos_byDL(int d, int line);
	int calcPos_byDR(int d, int row);

	int calcRow_byDL(int d, int line);
	int calcLine_byDR(int d, int row);
};

//-----------------------------------------------------------------------------
/** Матрица в диагональном формате. */
/** 0-я диагональ всегда главная диагональ. */
class MatrixDiagonal
{
public:
	typedef std::vector<real>::iterator iterator;
	typedef std::vector<real>::const_iterator const_iterator;

	//-------------------------------------------------------------------------
	MatrixDiagonal();
	MatrixDiagonal(int n, std::vector<int> format); // format[0] must be 0, because it's main diagonal
	MatrixDiagonal(const Matrix& a);

	void toDenseMatrix(Matrix& dense) const;
	void resize(int n, std::vector<int> format); // format[0] must be 0, because it's main diagonal

	void save(std::ostream& out) const;
	void load(std::istream& in);

	//-------------------------------------------------------------------------
	int dimension(void) const;
	int getDiagonalsCount(void) const;
	int getDiagonalSize(int diagNo) const;
	int getDiagonalPos(int diagNo) const;

	std::vector<int> getFormat(void) const;

	//-------------------------------------------------------------------------
	matrix_diagonal_iterator posBegin(int diagNo) const;
	matrix_diagonal_iterator posEnd(int diagNo) const;

	iterator begin(int diagNo);
	const_iterator begin(int diagNo) const;

	iterator end(int diagNo);
	const_iterator end(int diagNo) const;

private:
	std::vector<std::vector<real>> di;
	std::vector<int> fi;
	int n;
	Diagonal dc;
};

std::vector<int> makeSevenDiagonalFormat(int n, int m, int k);
std::vector<int> generateRandomFormat(int n, int diagonalsCount);

void generateDiagonalMatrix(
	int n,
	int min, int max,
	std::vector<int> format, 
	MatrixDiagonal& result
);

void generateDiagonallyDominantMatrix(
	int n, 
	std::vector<int> format, 
	bool isNegative,  
	MatrixDiagonal& result
);

bool mul(const MatrixDiagonal& a, const Vector& x, Vector& y);

//-----------------------------------------------------------------------------
/** Матричный "итератор" для движения по диагонали. */
class matrix_diagonal_iterator
{
public:
	matrix_diagonal_iterator(int n, int d, bool isEnd);

	matrix_diagonal_iterator& operator++();
	matrix_diagonal_iterator  operator++(int);

	bool operator==(const matrix_diagonal_iterator& b) const;
	bool operator!=(const matrix_diagonal_iterator& b) const;

	matrix_diagonal_iterator& operator+=(const ptrdiff_t& movement);

    int i, j;
};

/** Матричный "итератор" для движения по строке между различными диагоналями. Может обрабатывать как всю строку, так и только нижний треугольник. */
class matrix_diagonal_line_iterator
{
public:
	matrix_diagonal_line_iterator(int n, std::vector<int> format, bool isOnlyLowTriangle);

	matrix_diagonal_line_iterator& operator++();
	matrix_diagonal_line_iterator  operator++(int);

	bool isLineEnd(void) const;
	bool isEnd(void) const;

	int i, j; // i - текущая строка,  j - текущий столбец
	int d, dn, di; // d - формат текущей диагонали, d - номер текущей диагонали, di - номер текущего элемента в диагонали
private:
	std::map<int, int> m_map;
	std::vector<int> m_sorted_format;

	int line, start, end, pos;
	Diagonal dc;

	bool m_isLineEnd;
	bool m_isEnd;

	void calcPos(void);
};

//-----------------------------------------------------------------------------
/** Класс итеративного решателя СЛАУ для диагональной матрицы. */
struct IterationsResult
{
	int iterations;
	double relativeResidual;
};

class SolverSLAE_Iterative
{
public:
	SolverSLAE_Iterative();

	void save(std::ostream& out) const;
	void load(std::istream& in);

	IterationsResult jacobi(const MatrixDiagonal& a, const Vector& y, Vector& x) const;
	IterationsResult seidel(const MatrixDiagonal& a, const Vector& y, Vector& x) const;

	double 			w;
	bool 			isLog;
	std::ostream& 	log;
	Vector 			start;
	double 			epsilon;
	int 			maxIterations;

private:
	mutable Vector x1;

	// До итерации: x - текущее решение. После итерации x - следующее решение.
	void iteration_jacobi(const MatrixDiagonal& a, const Vector& y, Vector& x) const;
	void iteration_seidel(const MatrixDiagonal& a, const Vector& y, Vector& x) const;

	typedef std::function<void(const SolverSLAE_Iterative*, const MatrixDiagonal&, const Vector&, Vector&)> step_function;

	IterationsResult iteration_process(
		const MatrixDiagonal& a, 
		const Vector& y, 
		Vector& x, 
		step_function step
	) const;
};