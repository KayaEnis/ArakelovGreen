# Example 5.4

load('ArakelovGreen.sage')

var('p0','p1','p2','p3','a','b','c')
V = [p0, p1, p2, p3]
E = [(p0,p1,b), (p0,p2,a/2), (p2,p1,a/2), (p0,p3,c/2), (p3,p1,c/2)]
Γ = MetrizedGraph(V,E)

D = Γ.can_div()
true_value = 1/6*(a+b+c + a*b*c/(a*b+a*c+b*c))
our_value = Γ.eps_inv(D)
diff = (true_value - our_value).simplify_full()

print('------------------------------------------------------------------------------------')
print('The true value is')
print(true_value)
print()
print('Our value is')
print(our_value)
print()
print('The difference is ' + str(diff) + '. Success!')
print('------------------------------------------------------------------------------------')