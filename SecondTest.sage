# A simple script that computes the value matrix of a metrized graph and checks 
# if the values on vertices match the values obtained by direct computation. 

load('ArakelovGreen.sage')

# The Arakelov-Green function on vertices
def gmuD(Γ, D, p, q):
    n = Γ.order()

    tau = Γ.tau()

    gD = sum(D[i]*Γ.vol(p,q,Γ.v(i)) for i in range(n)) + 4*tau - Γ.res(p,q)
    gD *= 1/(sum(D)+2)

    cD = 8*tau*(sum(D)+1) + sum(D[i]*D[j]*Γ.res(Γ.v(i),Γ.v(j)) for i in range(n) for j in range(n))
    cD *= 1/(2*(sum(D)+2)^2)

    return gD - cD

# A big weird metrized graph with bridges, non-bridges, and different directions, 
# for illustration.
var('v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','w1','w2','w3','w4','w5','w6','w7','w8')
V = [v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, w1, w2, w3, w4, w5, w6, w7, w8]
E = [(v1,v2,12), (v1,v6,312), (v2,v7,123), (v4,v3,1321), (v4,v5,7777), (v4,v6,1), (v7,v6,11), (v7,v5,1321321), (v7,v8,2), (v8,v2,19), (v2,v3,8), (v8,v3,9), (v8,v4,1233), (v1,v9,2123), (v9,v10,321), (v10,w1,432), (w1,w2,1), (w2,w3,1), (w2,w4,2), (w3,w4,1), (w4,w1,1), (w2,w5,2), (w5,w6,10), (w6,w7,1), (w6,w8,4), (w7,w8,1)]
Γ = MetrizedGraph(V, E)
n = Γ.order()
m = Γ.size()

# Just some weird divisor, for the sake of illustration.
D = [(i^3 + 14) % 27 for i in range(n)]
ZD = Γ.val_mat(D)

test = []

for i in range(n):
    for j in range(n):
        v = Γ.v(i)
        w = Γ.v(j)

        direct_value = gmuD(Γ,D,v,w)

        F = []

        for nv in range(m):
            for nw in range(m):
                f = ZD[nv,nw]

                if v is Γ.p(nv) and w is Γ.p(nw):
                    F.append( f(x=0, y=0) )

                elif v is Γ.p(nv) and w is Γ.q(nw):
                    F.append( f(x=0, y=Γ.lth(nw)) )

                elif v is Γ.q(nv) and w is Γ.p(nw):
                    F.append( f(x=Γ.lth(nv), y=0) )

                elif v is Γ.q(nv) and w is Γ.q(nw):
                    F.append( f(x=Γ.lth(nv), y=Γ.lth(nw)) )

        if Set(F).cardinality() == 1 and F[0] == direct_value:
            test.append(0)
        else:
            print(F)
            test.append(1)

if sum(test) == 0:
    print('Success!')
else:
    print(test)
