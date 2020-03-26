from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import copy
import math

def split_numbers(sum_elem):
	number = parse_expr("1")
	other = parse_expr("1")
	for mul_element in sum_elem.args:
		if mul_element.is_Rational:
			number *= mul_element
		else:
			other *= mul_element
	return (number, other)

def insert_into_matrix_dict(M, key, value):
	if key not in M:
		M.update([(key, zeros(len(f)))])
	
	matrix = M.get(key)
	matrix[i, j] = matrix[i, j] + value
	M.update([(key, matrix)])

def lcm(a, b):
    return abs(a*b) // math.gcd(a, b)

def print_matrices(M, mul_to):
	for key, value in M.items():
		elems_lcm = 1
		for elem in value:
			elems_lcm = lcm(elems_lcm, fraction(elem)[1])
		key = key / elems_lcm
		value = simplify(value * elems_lcm)
		print("$${}{}\\cdot {} +$$".format(latex(mul_to), latex(key), latex(value)))

def get_function(list, var, xy, h):
	a, b, c = symbols("a b c")
	x, y = symbols("x y")
	a1 = solve([
		a*(xy +   0)*(xy +   0) + b*(xy +   0) + c - list[0], 
		a*(xy + h/2)*(xy + h/2) + b*(xy + h/2) + c - list[1], 
		a*(xy +   h)*(xy +   h) + b*(xy +   h) + c - list[2]
	], [a, b, c])
	return a1[a]*var*var + a1[b]*var + a1[c]

x, y = symbols("x y")
x_p, h_x = symbols("x_p h_x")
y_p, h_y = symbols("y_p h_y")

X1 = get_function([1, 0, 0], x, x_p, h_x)
X2 = get_function([0, 1, 0], x, x_p, h_x)
X3 = get_function([0, 0, 1], x, x_p, h_x)

Y1 = get_function([1, 0, 0], y, y_p, h_y)
Y2 = get_function([0, 1, 0], y, y_p, h_y)
Y3 = get_function([0, 0, 1], y, y_p, h_y)

print("X1: {}\nX2: {}\nX3: {}".format(X1, X2, X3))
print("Y1: {}\nY2: {}\nY3: {}".format(Y1, Y2, Y3))
print()

print("latex:")
print("X1: {}\nX2: {}\nX3: {}".format(latex(X1), latex(X2), latex(X3)))
print("Y1: {}\nY2: {}\nY3: {}".format(latex(Y1), latex(Y2), latex(Y3)))
print()

f = []
for i in [Y1, Y2, Y3]:
	for j in [X1, X2, X3]:
		f.append(j*i)

M = dict()
for i, a in enumerate(f):
	for j, b in enumerate(f):
		temp = simplify(
			integrate(
				integrate(
					a * b,
					(y, y_p, y_p + h_y)
				),
				(x, x_p, x_p + h_x)
			)
		)
		print("{}%".format((i * len(f) + j)/81.0 * 100.0))
		value, key = split_numbers(temp)
		insert_into_matrix_dict(M, key, value)

print("$$M = $$")
print_matrices(M, Symbol("gamma"))
print()

G = dict()
for i, a in enumerate(f):
	for j, b in enumerate(f):
		temp = simplify(
			integrate(
				integrate(
					diff(a, x) * diff(b, x) + diff(a, y) * diff(b, y),
					(y, y_p, y_p + h_y)
				),
				(x, x_p, x_p + h_x)
			)
		)
		print("{}%".format((i + j)/81.0 * 100.0))
		for sum_elem in temp.args:
			value, key = split_numbers(sum_elem)
			insert_into_matrix_dict(G, key, value)

print("$$G = $$")
print_matrices(G, Symbol("lambda"))
print()