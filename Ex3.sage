# Example 5.3

load('ArakelovGreen.sage')

# Initialize the tesseract graph
T = graphs.CubeGraph(4)

# Set edge lengths
for p,q,_ in T.edges():
    T.set_edge_label(p,q,1)

V = T.vertices()
E = T.edges()
Γ = MetrizedGraph(V, E)

# A silly divisor
D = [i for i in range(Γ.order())]
value_thm = Γ.eps_inv(D, method='thm')
value_lem = Γ.eps_inv(D, method='lem')

print('------------------------------------------------------------------------------------')
print('The value using "thm" is ' + str(value_thm))
print()
print('The value using "lem" is ' + str(value_lem))
print()
if value_thm == value_lem:
    print('Both methods yield the same value. Success!')
else:
    print('Values do not match, something is wrong!')
print('------------------------------------------------------------------------------------')
