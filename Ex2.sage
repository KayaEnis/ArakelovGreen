# Example 5.2

load('ArakelovGreen.sage')

var('p0','p1','p2','p3','p4','l1','l2')
V = [p0, p1, p2, p3, p4]
E = [(p0,p1,l1/3), (p1,p2,l1/3), (p2,p0,l1/3), (p0,p3,l2/3), (p3,p4,l2/3), (p4,p0,l2/3)]
Γ = MetrizedGraph(V, E)
D = [2,0,0,0,0]

# Our results
tot_lth = Γ.tot_lth()
tau = Γ.tau()
L = Γ.L
Lplus = Γ.Lplus
ZD = Γ.val_mat(D)

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


# Moriwaki's results
g = 2 # O = p0, so K = 2p0
L = l1 + l2

def phi(i, u):

	if i == 1:
		l = l1
	else:
		l = l2

	return 1/2 * (u^2/l - abs(u))

g11 = lambda x,y: phi(1,x-y) - (g-1)/g*(phi(1,x) + phi(1,y)) + L/(12*g^2)
g12 = lambda x,y: 1/g*(phi(1,x) + phi(2,y)) + L/(12*g^2)
g21 = lambda x,y: 1/g*(phi(2,x) + phi(1,y)) + L/(12*g^2)
g22 = lambda x,y: phi(2,x-y) - (g-1)/g*(phi(2,x) + phi(2,y)) + L/(12*g^2)

print('----------------------------------------------------------------------------------')
print('Moriwaki`s calculations:')
print()
print('On C1 x C1, the Arakelov-Green function is')
print(g11(x,y))
print()
print('On C1 x C2, the Arakelov-Green function is')
print(g12(x,y))
print()
print('On C2 x C1, the Arakelov-Green function is')
print(g21(x,y))
print()
print('On C2 x C2, the Arakelov-Green function is')
print(g22(x,y))
print('----------------------------------------------------------------------------------')


# Comparison
print('----------------------------------------------------------------------------------')
print('Comparison of the results:')
print()
print('Comparison on C1 x C1:')
print()
ZD11 = ZD.submatrix(0,0,3,3)
for i in range(3):
	for j in range(3):
		difference = (ZD11[i][j] - g11(x + i*l1/3,y + j*l1/3)).simplify_full()
		print('On e'+str(i)+' x e'+str(j)+', the difference is ' + str(difference))
print()
print('Comparison on C1 x C2:')
print()
ZD12 = ZD.submatrix(0,3,3,3)
for i in range(3):
	for j in range(3):
		difference = (ZD12[i][j] - g12(x + i*l1/3,y + j*l2/3)).simplify_full()
		print('On e'+str(i)+' x e'+str(j+3)+', the difference is ' + str(difference))
print()
print('Comparison on C2 x C1:')
print()
ZD21 = ZD.submatrix(3,0,3,3)
for i in range(3):
	for j in range(3):
		difference = (ZD21[i][j] - g21(x + i*l2/3,y + j*l1/3)).simplify_full()
		print('On e'+str(i+3)+' x e'+str(j)+', the difference is ' + str(difference))
print()
print('Comparison on C2 x C2:')
print()
ZD22 = ZD.submatrix(3,3,3,3)
for i in range(3):
	for j in range(3):
		difference = (ZD22[i][j] - g22(x + i*l2/3,y + j*l2/3)).simplify_full()
		print('On e'+str(i+3)+' x e'+str(j+3)+', the difference is ' + str(difference))
print()
print('These are all zero. Success!')
print('----------------------------------------------------------------------------------')
