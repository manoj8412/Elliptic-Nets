attach("path/EllNet.sage")
def NetSign(V, P, Q):

    if Q.curve() != P.curve():
	print("P and Q are from different curve...")
        return

    C = ComplexField(100)

    m = V[0]
    n = V[1]
    if m == 0 and n == 0:
        return 0
    E = P.curve()

    L = E.period_lattice()
    B1 = L.basis()[0]
    B2 = L.basis()[1]

    if B2.real() == 0:
        W1 = B1
        W2 = -B2

    else:
        B1 = 2*B2.real()
        B2 = B2
        W1 = B2
        W2 = -2*B2 + B1

    t = W1/W2
    q = e^(2*pi*I*t)


    u1 = exp(2*pi*I*P.elliptic_logarithm()/W2)
    while C(u1).real() < -1 or C(u1).real() > 1:
        u1 = u1*q

    u2 = exp(2*pi*I*Q.elliptic_logarithm()/W2)
    while C(u2).real() < -1 or C(u2).real() > 1:
	    u2 = u2*q

    b1 = log(abs(C(u1)), abs(C(q)))
    b2 = log(abs(C(u2)), abs(C(q)))

    if C(q).real() > 0:

        #print "q = ", C(q)
        #print "u1 = ", C(u1)
        #print "u2 = ", C(u2)


        if C(u1).real() > 0 and C(u2).real() > 0:
	        return (-1)^(Mod(floor(m*b1 + n*b2) + floor(b1 + b2)*m*n +m^2 +n^2 + m*n + 1,2))
        if C(u1).real() < 0 and C(u2).real() < 0:
	        if (m+n) % 2 == 0:
	            return (-1)^(Mod(floor(m*b1 + n*b2) + floor(b1 + b2)*m*n + floor(m/2) + floor(n/2) + m^2 + n^2 + m*n + 1,2))
	        else:
	            return (-1)^(Mod(floor(b1 + b2)*m*n +floor(m/2) + floor(n/2) + m^2 +n^2 + m*n + 1,2))
        if (C(u1)).real() < 0 and C(u2).real() > 0:
	        if m%2 == 0:
		        return (-1)^(Mod(floor(m*b1 + n*b2) + floor(b1 + b2)*m*n + m/2 + m^2 + n^2 + m*n + 1,2))
	        else:
		        return (-1)^(Mod(floor(b1 + b2)*m*n + (m-1)/2 +m^2 +n^2 + m*n + 1, 2))

        if (C(u1)).real() > 0 and C(u2).real() < 0:
	        if n%2 == 0:
		        return (-1)^(Mod(floor(m*b1 + n*b2) + floor(b1 + b2)*m*n + n/2 + m^2 + n^2 + m*n + 1,2))
	        else:
		        return (-1)^(Mod(floor(b1 + b2)*m*n + (n-1)/2 + m^2 +n^2 + m*n + 1, 2))

    if C(q).real() < 0:
	    while C(u1).real() < 0:
	       u1 = u1*q
	    while C(u2).real() < 0:
	       u2 = u2*q


	    #print "q = ", C(q)
	    #print "u1 = ", C(u1)
	    #print "u2 = ", C(u2)

	    b1 = log(abs(C(u1)),abs(C(q)))/2
	    b2 = log(abs(C(u2)),abs(C(q)))/2
	    return (-1)^(Mod(floor(m*b1+n*b2) + m^2 + n^2 + m*n + 1 + floor(b1+b2)*m*n, 2))
