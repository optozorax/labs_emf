#include "sparse.h"

//-----------------------------------------------------------------------------
matrix_sparse_t::matrix_sparse_t(int n) : n(n) {
	d.resize(n);
	i.resize(n+1, 0);
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::to_dense(matrix_t& m) const {
	m = matrix_t(n, n);
	m.fill(0);
	for (int _i = 0; _i < n; ++_i) {
		m(_i, _i) = d[_i];
		for (int _j = 0; _j < line_elem_count(_i); ++_j) {
			m(_i, line_elem_row(_i, _j)) = l[line_elem_start(_i) + _j];
			m(line_elem_row(_i, _j), _i) = u[line_elem_start(_i) + _j];
		}
	}
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::clear_line(int line) {
	d[line] = 0;
	for (int _i = 0; _i < i.size()-1; _i++) {
		for (int pj = i[_i]; pj < i[_i+1]; pj++) {
			int _j = j[pj];
			if (_j == line) u[pj] = 0;
			if (_i == line) l[pj] = 0;
		}
	}
}

//-----------------------------------------------------------------------------
int matrix_sparse_t::line_elem_start(int line) const { 
	return i[line]; 
}
//-----------------------------------------------------------------------------
int matrix_sparse_t::line_elem_row(int line, int elem) const { 
	return j[line_elem_start(line) + elem]; 
}

//-----------------------------------------------------------------------------
int matrix_sparse_t::line_elem_count(int line) const { 
	return i[line+1]-i[line]; 
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::decompose_lu_partial(matrix_sparse_t& lu) const {
	const matrix_sparse_t& a = *this;

	lu = a;
	for (int i = 0; i < lu.n; ++i) {
		// Заполняем нижний треугольник
		int line_start = lu.line_elem_start(i);
		int line_end = lu.line_elem_start(i+1);
		for (int j = line_start; j < line_end; ++j) {
			double sum = 0;

			int row = lu.j[j];
			int row_start = lu.line_elem_start(row);
			int row_end = lu.line_elem_start(row+1);

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
		int row_start = lu.line_elem_start(i);
		int row_end = lu.line_elem_start(i+1);
		for (int j = line_start; j < line_end; ++j) {
			double sum = 0;
			
			int line = lu.j[j];
			int line_start = lu.line_elem_start(line);
			int line_end = lu.line_elem_start(line+1);

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
		int line_row_start = lu.line_elem_start(i);
		int line_row_end = lu.line_elem_start(i+1);
		for (int j = line_row_start; j < line_row_end; ++j)
			sum += lu.l[j] * lu.u[j];

		lu.d[i] = sqrt(a.d[i] - sum);
	}
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::mul(vector_t& x_y) const {
	const matrix_sparse_t& a = *this;

	vector_t result(a.n); result.fill(0);

	for (int i = 0; i < a.n; ++i) {
		int start = a.line_elem_start(i);
		int size = a.line_elem_count(i);
		for (int j = 0; j < size; j++) {
			result[i] += a.l[start + j] * x_y[a.line_elem_row(i, j)];
			result[a.line_elem_row(i, j)] += a.u[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < a.n; ++i)
		result[i] += a.d[i] * x_y[i];

	x_y = result;
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::mul_t(vector_t& x_y) const {
	const matrix_sparse_t& a = *this;

	vector_t result(a.n); result.fill(0);

	for (int i = 0; i < a.n; ++i) {
		int start = a.line_elem_start(i);
		int size = a.line_elem_count(i);
		for (int j = 0; j < size; j++) {
			result(i) += a.u[start + j] * x_y[a.line_elem_row(i, j)];
			result(a.line_elem_row(i, j)) += a.l[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < a.n; ++i)
		result(i) += a.d[i] * x_y[i];

	x_y = result;
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::mul_l_inv_t(vector_t& x_y) const {
	const matrix_sparse_t& l = *this;

	for (int i = l.n - 1; i >= 0; i--) {
		int start = l.line_elem_start(i);
		int size = l.line_elem_count(i);

		x_y[i] /= l.d[i];
		for (int j = 0; j < size; ++j)
			x_y[l.line_elem_row(i, j)] -= x_y[i] * l.l[start + j];
	}
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::mul_u_inv_t(vector_t& x_y) const {
	const matrix_sparse_t& u = *this;
	
	for (int i = 0; i < u.n; ++i) {
		int start = u.line_elem_start(i);
		int size = u.line_elem_count(i);

		double sum = 0;
		for (int j = 0; j < size; ++j)
			sum += u.u[start + j] * x_y[u.line_elem_row(i, j)];
		x_y[i] = (x_y[i] - sum) / u.d[i];
	}
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::mul_l_inv(vector_t& x_y) const {
	const matrix_sparse_t& l = *this;
	
	for (int i = 0; i < l.n; ++i) {
		int start = l.line_elem_start(i);
		int size = l.line_elem_count(i);

		double sum = 0;
		for (int j = 0; j < size; ++j)
			sum += l.l[start + j] * x_y[l.line_elem_row(i, j)];
		x_y[i] = (x_y[i] - sum) / l.d[i];
	}
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::mul_u_inv(vector_t& x_y) const {
	const matrix_sparse_t& u = *this;
	
	for (int i = u.n-1; i >= 0; i--) {
		int start = u.line_elem_start(i);
		int size = u.line_elem_count(i);

		x_y[i] /= u.d[i];
		for (int j = 0; j < size; ++j)
			x_y[u.line_elem_row(i, j)] -= x_y[i] * u.u[start + j];
	}
}

//-----------------------------------------------------------------------------
void matrix_sparse_t::mul_u(vector_t& x_y) const {
	const matrix_sparse_t& u = *this;
	
	vector_t result(u.n); result.fill(0);

	for (int i = 0; i < u.n; ++i) {
		int start = u.line_elem_start(i);
		int size = u.line_elem_count(i);
		for (int j = 0; j < size; j++) {
			result[u.line_elem_row(i, j)] += u.u[start + j] * x_y[i];
		}
	}

	// Умножение диагональных элементов на вектор
	for (int i = 0; i < u.n; ++i)
		result[i] += u.d[i] * x_y[i];

	x_y = result;
}

//-----------------------------------------------------------------------------
ostream& operator<<(ostream& out, const matrix_sparse_t& m) {
	matrix_t dense;
	m.to_dense(dense);
	out << dense;
	return out;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
matrix_sparse_ra_t::matrix_sparse_ra_t(int n) : n(n), dm(n, 0), lm(n), um(n) {}

//-----------------------------------------------------------------------------
double& matrix_sparse_ra_t::operator()(int i, int j) {
	if (i == j) {
		return dm[i];
	} else if (i > j) {
		um[i][j] += 0;
		return lm[i][j];
	} else {
		lm[j][i] += 0;
		return um[j][i];
	}
}

//-----------------------------------------------------------------------------
const double& matrix_sparse_ra_t::operator()(int i, int j) const {
	if (i == j) {
		return dm[i];
	} else if (i > j) {
		if (lm[i].find(j) != lm[i].end())
			return lm[i].at(j);
	} else {
		if (um[j].find(i) != um[i].end())
			return um[j].at(j);
	}

	throw exception();
}

//-----------------------------------------------------------------------------
matrix_sparse_t matrix_sparse_ra_t::to_sparse(void) const {
	matrix_sparse_t result(n);
	result.n = dm.size();
	result.d = dm;
	for (int i = 0; i < lm.size(); ++i) {
		result.i[i+1] = result.i[i] + lm[i].size();
		for (auto& j : lm[i]) {
			result.j.push_back(j.first);
			result.l.push_back(j.second);
			result.u.push_back(um[i].at(j.first));
		}
	}
	return result;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
void mul(const vector_t& d, vector_t& x_y) {
	for (int i = 0; i < d.size(); i++)
		x_y[i] *= d[i];
}

//-----------------------------------------------------------------------------
void mul_inv(const vector_t& d, vector_t& x_y) {
	for (int i = 0; i < d.size(); i++)
		x_y[i] /= d[i];
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
vector_t solve_by_los_lu(
	const matrix_sparse_t& a,
	const vector_t& b,
	int maxiter,
	double eps,
	bool is_log
) {
	matrix_sparse_t lu(a.n);
	vector_t r, z, p;
	vector_t x, t1, t2;

	int n = a.n;

	a.decompose_lu_partial(lu);
	x = vector_t(n); x.fill(0);

	r = x;
	a.mul(r);
	for (int i = 0; i < n; i++)
		r[i] = b[i] - r[i];
	lu.mul_l_inv(r);

	z = r;
	lu.mul_u_inv(z);

	p = z;
	a.mul(p);
	lu.mul_l_inv(p);

	double flen = sqrt(b.dot(b));
	double residual;

	int i = 0;
	while (true) {
		double pp = p.dot(p);
		double alpha = (p.dot(r)) / pp;
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * z[i];
			r[i] -= alpha * p[i];
		}
		t1 = r;
		lu.mul_u_inv(t1);
		t2 = t1;
		a.mul(t2);
		lu.mul_l_inv(t2);
		double beta = -(p.dot(t2)) / pp;
		for (int i = 0; i < n; ++i) {
			z[i] = t1[i] + beta * z[i];
			p[i] = t2[i] + beta * p[i];
		}
		residual = r.norm() / flen;
		i++;

		//if (is_log) cout << "Iteration: " << setw(4) << i << ", Residual: " << setw(20) << setprecision(16) << residual << endl;
		if (fabs(residual) < eps || i > maxiter)
			break;
	}

	return x;
}