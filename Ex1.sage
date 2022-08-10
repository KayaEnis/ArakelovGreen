# Example 5.1

load('ArakelovGreen.sage')

var('p0','p1','p2')
a = 1/2*pi
b = 1/2*pi
c = pi
V = [p1 ,p0, p2]
E = [(p0,p1,a), (p0,p2,c), (p1,p2,b)]
Γ = MetrizedGraph(V,E)

# Our results
tot_lth = Γ.tot_lth()
tau = Γ.tau()
L = Γ.L
Lplus = Γ.Lplus
ZD = Γ.val_mat()

print('------------------------------------------------------------------------------------')
print('Our calculations:')
print()
print('The total length is ' + str(tot_lth))
print()
print('The tau constant is ' + str(tau))
print()
print('The discrete Laplacian matrix is')
print(L)
print()
print('The Moore--Penrose generalized inverse is')
print(Lplus)
print()
print('The value matrix is')
print(ZD)
print('------------------------------------------------------------------------------------')

# Cinkir's results
tot_lthC = a + b + c
tauC = tot_lthC/12
LC = Matrix([[1/a+1/b, -1/a, -1/b],[-1/a, 1/a+1/c, -1/c],[-1/b, -1/c, 1/b+1/c]])
LplusC = 1/(9*tot_lthC) * Matrix([[b*c + a*(4*b+c), b*c - 2*a*(b+c), -2*b*c + a*(-2*b+c)],[b*c - 2*a*(b+c), b*c + a*(b+4*c), a*(b-2*c) - 2*b*c],[-2*b*c + a*(-2*b+c), a*(b-2*c) - 2*b*c, 4*b*c + a*(b+c)]])
ZDC = tau*matrix.ones(3) - 1/(2*tot_lthC) * Matrix([[-(x-y)^2 + (a+b+c)*abs(x-y), (a+b+c-x-y)*(x+y), (b+c+x-y)*(a-x+y)],[(a+b+c-x-y)*(x+y), -(x-y)^2 + (a+b+c)*abs(x-y), (b+c-x-y)*(a+x+y)],[(b+c+x-y)*(a-x+y), (b+c-x-y)*(a+x+y), -(x-y)^2 + (a+b+c)*abs(x-y)]])

print('------------------------------------------------------------------------------------')
print('Cinkir`s calculations:')
print()
print('The total length is ' + str(tot_lthC))
print()
print('The tau constant is ' + str(tauC))
print()
print('The discrete Laplacian matrix is')
print(LC)
print()
print('The Moore--Penrose generalized inverse is')
print(LplusC)
print()
print('The value matrix is')
print(ZDC)
print('------------------------------------------------------------------------------------')


# Comparison
print('------------------------------------------------------------------------------------')
print('Comparison of the results:')
print()
print('The difference of the total lengths is ' + str(tot_lth - tot_lthC))
print()
print('The difference of the tau constants is ' + str(tau - tauC))
print()
print('The difference of the discrete Laplacian matrices is')
print(L - LC)
print()
print('The difference of the Moore--Penrose generalized inverses is')
print(Lplus - LplusC)
print()
print('The difference of the value matrices is')
print((ZD - ZDC).simplify_full())
print('------------------------------------------------------------------------------------')

# Checking symmetry
ZDt = ZD.transpose().substitute(x=y,y=x)
ZDCt = ZDC.transpose().substitute(x=y,y=x)

print('------------------------------------------------------------------------------------')
print('Checking symmetry:')
print()
if (ZD - ZDt).simplify_full() == 0:
	print('Our matrix is symmetric. Success!')
else:
	print('Our matrix is not symmetric. Fail!')
print()
if (ZDC - ZDCt).simplify_full() == 0:
	print('Cinkir`s matrix is symmetric. Success!')
else:
	print('Cinkir`s matrix is not symmetric. Fail!')
print('------------------------------------------------------------------------------------')
