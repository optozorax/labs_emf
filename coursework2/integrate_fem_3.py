from sympy import *
import copy

a, b, c = symbols("a b c")

x, y = symbols("x y")

def get_function(list, var, xy, h):
	a1 = solve([
		a*xy*xy + b*xy + c - list[0], 
		a*(xy + h/2)*(xy + h/2) + b*(xy + h/2) + c - list[1], 
		a*(xy + h)*(xy + h) + b*(xy + h) + c - list[2]
	], (a, b, c))
	return a1[a]*var*var + a1[b]*var + a1[c]

x_p, h_x = symbols("x_p h_x")
x_p1 = x_p + h_x

y_p, h_y = symbols("y_p h_y")
y_p1 = y_p + h_y

X1 = get_function([1, 0, 0], x, x_p, h_x)
X2 = get_function([0, 1, 0], x, x_p, h_x)
X3 = get_function([0, 0, 1], x, x_p, h_x)

Y1 = get_function([1, 0, 0], y, y_p, h_y)
Y2 = get_function([0, 1, 0], y, y_p, h_y)
Y3 = get_function([0, 0, 1], y, y_p, h_y)

f = []

for i in [Y1, Y2, Y3]:
	for j in [X1, X2, X3]:
		f.append(j*i)

G1 = zeros(len(f))
G2 = zeros(len(f))

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
		(first, second) = temp.args
		(coef, up, down) = first.args
		if up == h_x:
			G1[i, j], G2[i, j] = first, second
		else:
			G1[i, j], G2[i, j] = second, first

M = zeros(len(f))
for i, a in enumerate(f):
	for j, b in enumerate(f):
		M[i, j] = simplify(
			integrate(
				integrate(
					a * b,
					(y, y_p, y_p + h_y)
				),
				(x, x_p, x_p + h_x)
			)
		)

first = h_x/(h_y*6*15)
second = h_y/(h_x*6*15)
G1 = simplify(G1 / first)
G2 = simplify(G2 / second)
l = Symbol("lambda")

print("G = {}\\cdot{} + {}\\cdot{}".format(latex(l * first), latex(G1), latex(l * second), latex(G2)))

multiplier = h_x*h_y/36/25
M = simplify(M / multiplier)
g = Symbol("gamma")

print("M = {}\\cdot{}".format(latex(g * multiplier), latex(M)))
