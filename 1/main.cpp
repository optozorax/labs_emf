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

	bool isBorder;
	//int conditionType; // 0 - первое, 1 - второе, 2 - третье

	double value;
	double x, y;

	Cell *up, *down, *left, *right;
};

typedef std::vector<Cell> Cells;

class Field
{
public:
	virtual bool isPointInside(double x, double y) const = 0;
};

class ShField : public Field
{
public:
	ShField(double sizex, double sizey) : sizex(sizex), sizey(sizey) {
	}
	bool isPointInside(double x, double y) const {
		// TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
		if (x >= 0 && x <= sizex && y >= 0 && y <= sizey)
			return true;
		else
			return true;
	}
private:
	double sizex, sizey;
};

class Grid
{
public:
	virtual Cells makeCells(const Field& field, 
	                        double startx, double starty,
	                        double sizex, double sizey,
	                        int isizex, int isizey) = 0;
};

class UniformGrid : public Grid
{
public:
	Cells makeCells(const Field& field, 
	                double startx, double starty,
	                double sizex, double sizey,
	                int isizex, int isizey) {
		Cells cells;
		cells.reserve(isizex * isizey);
		// [строка][столбец]
		std::vector<std::vector<int>> grid(isizey, std::vector<int>(isizex, -1));
		double x = startx, y = starty;
		double hx = sizex/(isizex-1), hy = sizey/(isizey-1);
		int k = 0;
		for (int i = 0; i < isizey; ++i, y += hy) {
			x = startx;
			for (int j = 0; j < isizex; ++j, x += hx) {		
				if (field.isPointInside(x, y)) {
					grid[i][j] = k;
					cells.push_back({k, false, 0, x, y,
						nullptr,
						(i != 0) ? &cells[grid[i-1][j]] : nullptr, 
						(j != 0) ? &cells[grid[i][j-1]] : nullptr,
						nullptr,
					});
					k++;
					if (i != 0) cells[grid[i-1][j]].up = &cells.back();
					if (j != 0) cells[grid[i][j-1]].right = &cells.back();
				}
			}
		}

		return cells;
	}
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
					cells.push_back({k, false, 0, x, y,
						nullptr,
						(i != 0) ? &cells[grid[i-1][j]] : nullptr, 
						(j != 0) ? &cells[grid[i][j-1]] : nullptr,
						nullptr,
					});
					k++;
					if (i != 0) cells[grid[i-1][j]].up = &cells.back();
					if (j != 0) cells[grid[i][j-1]].right = &cells.back();
				}
			}
		}

		return cells;
	}
private:
	double c;
};

void fillWithFunction(Cells& cells, const Function2D& f) {
	for (auto& i : cells)
		i.value = f(i.x, i.y);
}

void fillBoundaryConditions1(Cells& cells, const Function2D& f) {
	// Первые краевые условия
	for (auto& i : cells) {
		if (i.up == nullptr || i.down == nullptr || i.left == nullptr || i.right == nullptr) {
			i.isBorder = true;
			i.value = f(i.x, i.y);
		}
	}
}

void fillBoundaryConditions3(Cells& cells, const Function2D& f) {

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
		if (cell.isBorder) {
			pushLine(A, i, {{cell.i, 1}});
			b[i] = cell.value;
		} else
		/*if (cell.up == nullptr) {
			
		} else if (cell.down == nullptr) {

		} else if (cell.left == nullptr) {

		} else if (cell.right == nullptr) {

		} else */
		{
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

#ifdef _DEBUG
#define debug(a) std::cout << #a << ": " << std::endl << (a) << std::endl;
#else
#define debug(a) ;
#endif

int main() {
	std::ofstream fout("a.txt");

	for (int i = 6; i <= 300; i+=2) {
		if (i % 10 == 0) std::cout << i << std::endl;
		ShField field(1, 1);
		//UniformGrid grid;
		NonUniformGrid grid(i/100.0);
		int size = 15;
		auto cells_answer = grid.makeCells(field, 0, 0, 1, 1, size, size);
		auto cells_question = grid.makeCells(field, 0, 0, 1, 1, size, size);
		auto f = [] (double x, double y) -> double { return exp(x*y + x*x*y + 3); };
		auto rightPart = [] (double x, double y) -> double { return exp(x*y + x*x*y + 3)*(2*x*x*x + x*x*x*x + 4*x*y*y + y*(2+y) + x*x*(1+4*y*y)); };
		/*auto f = [] (double x, double y) -> double { return exp(x*y); };
		auto rightPart = [] (double x, double y) -> double { return exp(x*y)*(x*x+y*y); };*/
		/*auto f = [] (double x, double y) -> double { return x*x*x*x + y*y*y*y; };
		auto rightPart = [] (double x, double y) -> double { return 12*(x*x+y*y); };*/
		/*auto f = [] (double x, double y) -> double { return x*x*x + y*y*y; };
		auto rightPart = [] (double x, double y) -> double { return 6*(x+y); };*/
		/*auto f = [] (double x, double y) -> double { return x*x + y*y; };
		auto rightPart = [] (double x, double y) -> double { return 4; };*/
		/*auto f = [] (double x, double y) -> double { return 2*x + y; };
		auto rightPart = [] (double x, double y) -> double { return 0; };*/
		fillWithFunction(cells_answer, f);
		fillBoundaryConditions1(cells_question, f);
		auto slae = makeSLAE(cells_question, rightPart);
		//debug(slae.first);
		//debug(slae.second);
		auto answer = solveSLAE(slae); //debug(answer);
		setCells(cells_question, answer);
		double difference = calcDifference(cells_answer, cells_question);
		//std::cout << "Precision of SLAE solve: " << (slae.first*answer - slae.second).norm() << std::endl;
		//std::cout << "Difference: " << difference << std::endl;
		fout << i/100.0 << "\t" << difference << std::endl;
	}

	fout.close();

	system("pause");
}