#include <iostream>
#include <functional>
#include <fstream>
#include <Eigen/Eigen>
#include <Eigen/Dense>

typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;

typedef std::function<double(double)> Function1D;
typedef std::function<double(double, double)> Function2D;

struct Cell
{
	int i;

	enum 
	{
		INNER, // Нет краевых условий
		FIRST, // f(x, y) = value
		SECOND,  // f'(x, y) = value
		THIRD,  // f'(x, y) + c1 * f(x, y) + c2 = 0
	}
	conditionType; // краевое условие

	union {
		struct {
			double value;
		} first;
		struct {
			Cell* cell_around;
			double value;
			double h;
		} second;
		struct {
			Cell* cell_around;
			double value;
			double h;
			double c1;
		} third;
	} conditionValue;

	double value;
	double x, y;

	Cell *up, *down, *left, *right;
};

typedef std::vector<Cell> Cells;

class Field
{
public:
	Field(double startx, double starty, double sizex, double sizey) : sizex(sizex), sizey(sizey), startx(startx), starty(starty) {
	}

	virtual bool isPointInside(double x, double y) const = 0;

	double sizex, sizey;
	double startx, starty;
};

bool isInsideSquare(double x, double y, double sizex, double sizey, double startx, double starty) {
	x -= startx;
	y -= starty;
	return x >= 0 && x <= sizex && y >= 0 && y <= sizey;
}

class SquareField : public Field
{
public:
	SquareField(double startx, double starty, double sizex, double sizey) : Field(startx, starty, sizex, sizey) {
	}
	bool isPointInside(double x, double y) const {
		if (isInsideSquare(x, y, sizex, sizey, startx, starty))
			return true;
		else
			return false;
	}
};

class ShField : public Field
{
public:
	ShField(double startx, double starty, double sizex, double sizey) : Field(startx, starty, sizex, sizey) {
	}
	bool isPointInside(double x, double y) const {
		if (isInsideSquare(x, y, sizex, sizey, startx, starty) && 
		    !isInsideSquare(x, y, sizex * 0.2, sizey * 0.2, startx + sizex * 0.2, starty + sizey * 0.8) &&
		    !isInsideSquare(x, y, sizex * 0.2, sizey * 0.2, startx + sizex * 0.6, starty + sizey * 0.8))
			return true;
		else
			return false;
	}
};

class Grid
{
public:
	virtual Cells makeCells(const Field& field, 
	                        double startx, double starty,
	                        double sizex, double sizey,
	                        int isizex, int isizey) = 0;
};

class NonUniformGrid : public Grid 
{
public:
	NonUniformGrid(double c) : c(c) {}

	Cells makeCells(const Field& field, 
	                double startx, double starty,
	                double sizex, double sizey,
	                int isizex, int isizey) {
		Cells cells;
		cells.reserve(isizex * isizey);
		// [строка][столбец]
		std::vector<std::vector<int>> grid(isizey, std::vector<int>(isizex, -1));

		double x = startx, y = starty;
		double hx, hy;
		if (std::fabs(c - 1) > 0.00001) {
			hx = sizex * (1.0-c)/(1-pow(c, isizex));
			hy = sizey * (1.0-c)/(1-pow(c, isizey));
		} else {
			// Это равномерная сетка, что ты здесь забыл?
			hx = sizex/(isizex-1); 
			hy = sizey/(isizey-1);
		}
		double oldhx = hx;
		int k = 0;
		for (int i = 0; i < isizey; ++i, y += hy, hy *= c) {
			x = startx;
			hx = oldhx;
			for (int j = 0; j < isizex; ++j, x += hx, hx *= c) {		
				if (field.isPointInside(x, y)) {
					grid[i][j] = k;
					cells.push_back({k, Cell::INNER, 0, 0, x, y,
									nullptr,
									(i != 0 && grid[i-1][j] != -1) ? &cells[grid[i-1][j]] : nullptr, 
									(j != 0 && grid[i][j-1] != -1) ? &cells[grid[i][j-1]] : nullptr,
									nullptr,
									});
					k++;
					if (i != 0 && grid[i-1][j] != -1) cells[grid[i-1][j]].up = &cells.back();
					if (j != 0 && grid[i][j-1] != -1) cells[grid[i][j-1]].right = &cells.back();
				}
			}
		}

		return cells;
	}
private:
	double c;
};

class UniformGrid : public Grid
{
public:
	UniformGrid() : grid(1.0) {
	}

	Cells makeCells(const Field& field, 
	                double startx, double starty,
	                double sizex, double sizey,
	                int isizex, int isizey) {
		return grid.makeCells(field, startx, starty, sizex, sizey, isizex, isizey);
	}
private:
	NonUniformGrid grid;
};

Function2D calcFirstDerivativeX(const Function2D& f) {
	return [f] (double x, double y) -> double {
		const double h = 0.0001;
		return (-f(x+2*h, y)+8*f(x+h, y)-8*f(x-h, y)+f(x-2*h, y))/(12*h);
	};
}

Function2D calcFirstDerivativeY(const Function2D& f) {
	return [f] (double x, double y) -> double {
		const double h = 0.0001;
		return (-f(x, y+2*h)+8*f(x, y+h)-8*f(x, y-h)+f(x, y-2*h))/(12*h);
	};
}

Function2D calcSecondDerivativeX(const Function2D& f) {
	return [f] (double x, double y) -> double {
		const double h = 0.001;
		return (-f(x+2*h, y) + 16*f(x+h, y) - 30*f(x, y) + 16*f(x-h, y) - f(x-2*h, y))/(12*h*h);
	};
}

Function2D calcSecondDerivativeY(const Function2D& f) {
	return [f] (double x, double y) -> double {
		const double h = 0.001;
		return (-f(x, y+2*h) + 16*f(x, y+h) - 30*f(x, y) + 16*f(x, y-h) - f(x, y-2*h))/(12*h*h);
	};
}

Function2D calcLaplacian(const Function2D& f) {
	auto fxx = calcSecondDerivativeX(f);
	auto fyy = calcSecondDerivativeY(f);
	return [fxx, fyy] (double x, double y) -> double {
		return fxx(x, y) + fyy(x, y);
	};
}

void fillWithFunction(Cells& cells, const Function2D& f) {
	for (auto& i : cells) {
		i.value = f(i.x, i.y);
	}
}

void fillBoundaryConditions1(Cells& cells, const Function2D& f) {
	// Первые краевые условия
	for (auto& i : cells) {
		if (i.up == nullptr || i.down == nullptr || i.left == nullptr || i.right == nullptr) {
			i.conditionType = Cell::FIRST;
			i.conditionValue.first.value = f(i.x, i.y);
		}
	}
}

void fillBoundaryConditions2(Cells& cells, const Function2D& f) {
	auto fx = calcFirstDerivativeX(f);
	auto fy = calcFirstDerivativeY(f);
	int isFirst = 2;
	for (auto& i : cells) {
		if (i.up == nullptr || i.down == nullptr || i.left == nullptr || i.right == nullptr) {
			if (isFirst || i.left == nullptr && i.right == nullptr || i.up == nullptr && i.down == nullptr) {
				// Тут уже ничего не сделаешь, потому что это палка, у которой нет соседей, чтобы можно было сделать краевые условия с производной
				i.conditionType = Cell::FIRST;
				i.conditionValue.first.value = f(i.x, i.y);
				isFirst--;
			} else {
				i.conditionType = Cell::SECOND;
				double fxv = fx(i.x, i.y);
				double fyv = fy(i.x, i.y);
				
				Cell* cell_around;
				double& result = i.conditionValue.second.value;
				if (i.left == nullptr) {
					result = -fxv;
					cell_around = i.right;
				} else if (i.right == nullptr) {
					result = fxv;
					cell_around = i.left;
				} else if (i.up == nullptr) {
					result = fyv;
					cell_around = i.down;
				} else if (i.down == nullptr) {
					result = -fyv;
					cell_around = i.up;
				} else {
					throw std::exception("This is inner cell!");
				}
				i.conditionValue.second.cell_around = cell_around;

				double hx = std::fabs(cell_around->x - i.x);
				double hy = std::fabs(cell_around->y - i.y);
				double h = std::max(hx, hy);
				i.conditionValue.second.h = h;
			}
		}
	}
}

void fillBoundaryConditions3(Cells& cells, const Function2D& f, double c1) {
	auto fx = calcFirstDerivativeX(f);
	auto fy = calcFirstDerivativeY(f);
	for (auto& i : cells) {
		if (i.up == nullptr || i.down == nullptr || i.left == nullptr || i.right == nullptr) {
			if (i.left == nullptr && i.right == nullptr || i.up == nullptr && i.down == nullptr) {
				// Тут уже ничего не сделаешь, потому что это палка, у которой нет соседей, чтобы можно было сделать краевые условия с производной
				i.conditionType = Cell::FIRST;
				i.conditionValue.first.value = f(i.x, i.y);
			} else {
				i.conditionType = Cell::THIRD;
				double fxv = fx(i.x, i.y);
				double fyv = fy(i.x, i.y);

				Cell* cell_around;
				double& result = i.conditionValue.third.value;
				if (i.left == nullptr) {
					result = -fxv;
					cell_around = i.right;
				} else if (i.right == nullptr) {
					result = fxv;
					cell_around = i.left;
				} else if (i.up == nullptr) {
					result = fyv;
					cell_around = i.down;
				} else if (i.down == nullptr) {
					result = -fyv;
					cell_around = i.up;
				} else {
					throw std::exception("This is inner cell!");
				}
				result += f(i.x, i.y)*c1;
				i.conditionValue.third.cell_around = cell_around;

				double hx = std::fabs(cell_around->x - i.x);
				double hy = std::fabs(cell_around->y - i.y);
				double h = std::max(hx, hy);
				i.conditionValue.third.h = h;
				i.conditionValue.third.c1 = c1;
			}
		}
	}
}

double calcDifference(const Cells& a, const Cells& b) {
	if (a.size() != b.size())
		return -1;
	double sum = 0;
	for (int i = 0; i < a.size(); ++i) {
		if (a[i].x != b[i].x || a[i].y != b[i].y)
			return -1;
		sum += (a[i].value - b[i].value)*(a[i].value - b[i].value);
	}
	return std::sqrt(sum);
}

struct PushType
{
	int i;
	double value;
};

void pushLine(Matrix& A, int line, std::vector<PushType> p) {
	std::sort(p.begin(), p.end(), [] (auto& a, auto& b) -> bool {
		return a.i < b.i;
	});

	for (int i = 0; i < A.cols(); ++i) {
		for (auto& j : p) {
			if (i == j.i) {
				A(line, i) = j.value;
				goto next_line;
			}
		}
		A(line, i) = 0;

		next_line:;
	}
}

std::pair<Matrix, Vector> makeSLAE(const Cells& cells, const Function2D& rightPart) {
	Matrix A(cells.size(), cells.size());
	Vector b(cells.size());

	for (int i = 0; i < cells.size(); ++i) {
		const auto& cell = cells[i];
		switch (cell.conditionType) {
			case Cell::INNER: {
				double uph = std::fabs(cell.up->y - cell.y);
				double downh = std::fabs(cell.down->y - cell.y);
				double lefth = std::fabs(cell.left->x - cell.x);
				double righth = std::fabs(cell.right->x - cell.x);
				pushLine(A, i, {
					{cell.left->i, 2.0/(lefth*(righth + lefth))}, 
					{cell.down->i, 2.0/(downh*(uph + downh))}, 
					{cell.right->i, 2.0/(righth*(righth + lefth))}, 
					{cell.up->i, 2.0/(uph*(uph + downh))}, 
					{cell.i, -2.0/(uph*downh)-2.0/(lefth*righth)}, 
				});

				b[i] = rightPart(cell.x, cell.y);
			} break;
			case Cell::FIRST: {
				pushLine(A, i, {{cell.i, 1}});
				b[i] = cell.conditionValue.first.value;
			} break;
			case Cell::SECOND: {
				double h = cell.conditionValue.second.h;
				pushLine(A, i, {
					{cell.i, 1.0/h}, 
					{cell.conditionValue.second.cell_around->i, -1.0/h}, 
				});
				b[i] = cell.conditionValue.second.value;
			} break;
			case Cell::THIRD: {
				double h = cell.conditionValue.third.h;
				pushLine(A, i, {
					{cell.i, 1.0/h + cell.conditionValue.third.c1}, 
					{cell.conditionValue.third.cell_around->i, -1.0/h}, 
				});
				b[i] = cell.conditionValue.second.value;
			} break;
		}
	}

	return {A, b};
}

Vector solveSLAE(std::pair<Matrix, Vector>& slae) {
	Matrix& A = slae.first;
	Vector& b = slae.second;
	using namespace Eigen;
	return A.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
}

void setCells(Cells& cells, const Vector& answer) {
	for (int i = 0; i < cells.size(); i++)
		cells[i].value = answer(i);
}

void make_table_nonuniform_grid(const std::string& name, const Field& field, const Function2D& f, int size, int boundaryCondition) {
	auto rightPart = calcLaplacian(f);
	std::ofstream fout(name);
	fout << "c\tnorm" << std::endl;
	for (int i = 1; i <= 400; i++) {
		NonUniformGrid grid(i/100.0);
		auto cells_answer = grid.makeCells(field, field.startx, field.starty, field.sizex, field.sizey, size, size);
		auto cells_question = grid.makeCells(field, field.startx, field.starty, field.sizex, field.sizey, size, size);
		fillWithFunction(cells_answer, f);
		if (boundaryCondition == 1) {
			fillBoundaryConditions1(cells_question, f);
		} else if (boundaryCondition == 2) {
			fillBoundaryConditions2(cells_question, f);
		} else {
			fillBoundaryConditions3(cells_question, f, 1);
		}
		auto slae = makeSLAE(cells_question, rightPart);
		auto answer = solveSLAE(slae);
		setCells(cells_question, answer);
		double difference = calcDifference(cells_answer, cells_question);
		if (difference == difference)
			fout << i/100.0 << "\t" << difference << std::endl;

		std::cout << "\r" << i;
	}
	fout.close();
}

void make_table_size(const std::string& name, const Field& field, const Function2D& f, int boundaryCondition) {
	auto rightPart = calcLaplacian(f);
	std::ofstream fout(name);
	fout << "size\tnorm" << std::endl;
	for (int i = 5; i <= 50; i++) {
		UniformGrid grid;
		auto cells_answer = grid.makeCells(field, field.startx, field.starty, field.sizex, field.sizey, i, i);
		auto cells_question = grid.makeCells(field, field.startx, field.starty, field.sizex, field.sizey, i, i);
		fillWithFunction(cells_answer, f);
		if (boundaryCondition == 1) {
			fillBoundaryConditions1(cells_question, f);
		} else if (boundaryCondition == 2) {
			fillBoundaryConditions2(cells_question, f);
		} else {
			fillBoundaryConditions3(cells_question, f, 1);
		}
		auto slae = makeSLAE(cells_question, rightPart);
		auto answer = solveSLAE(slae);
		setCells(cells_question, answer);
		double difference = calcDifference(cells_answer, cells_question);
		if (difference == difference)
			fout << i << "\t" << difference << std::endl;

		std::cout << "\r" << i;
	}
	fout.close();
}

int main() {
	auto f = [] (double x, double y) -> double { return exp(x*y + x*x*y + 3); };
	//auto f = [] (double x, double y) -> double { return exp(x*y); };
	//auto f = [] (double x, double y) -> double { return sin(x) + cos(y); };
	//auto f = [] (double x, double y) -> double { return x*x*x*x*x + y*y*y*y*y; };
	//auto f = [] (double x, double y) -> double { return x*x*x*x + y*y*y*y; };
	//auto f = [] (double x, double y) -> double { return x*x*x + y*y*y; };
	//auto f = [] (double x, double y) -> double { return x*x + y*y; };
	//auto f = [] (double x, double y) -> double { return 2*x + y; };

	SquareField field(0, 0, 1, 1);
	//ShField field(0, 0, 1, 1);
	//make_table_size("a.txt", field, f);
	make_table_nonuniform_grid("a.txt", field, f, 14, 3);
}