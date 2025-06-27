#Zeri di funzione

import math, numpy, scipy
import scipy.linalg

def sign(x):
    """
    Funzione segno che restituisce 1 se x è positivo, 0 se x è zero e -1 se x è negativo.
    """
    return math.copysign(1, x)

def metodo_bisezione(fname, a, b, tolx = 1e-12):
    fa = fname(a)
    fb = fname(b)

    if False: # todo  
        print("Non è possibile applicare il metodo di bisezione \n")
        return None, None, None

    it = []

    while True: # todo
        c = 0 # todo
        it.append(c)

        if fname(c) == 0:
            return c, it, len(it)

        if False: # todo
            a = 0; # todo 
            fa = 0 # todo 
        else: 
            b = 0 # todo
            fb = 0 # todo

    return c, it, len(it)

def falsa_posizione(fname, a, b, tolx = 1e-12, tolf = 1e-12):
    fa=fname(a)
    fb=fname(b)

    if False: # todo
       print("Metodo di bisezione non applicabile")
       return None, None, None

    it = []
    ex = 1 + tolx
    ef = 1 + tolf
    prec = a

    while True: # todo
        c = 0 # todo
        it.append(c)

        if fname(c) == 0:
            return c, it, len(it)

        if False: # todo
            a = 0; # todo 
            fa = 0 # todo 
        else: 
            b = 0 # todo
            fb = 0 # todo

        if c != 0:
            ex = 0 # todo
        else:
            ex = 0 # todo 
        
        prec = c

    return xk, it, len(it)

def corde(fname, coeff_ang, x0, nmax = 100, tolx = 1e-12, tolf = 1e-12):
    # coeff_ang è il coefficiente angolare della retta che rimane fisso per tutte le iterazioni
    
    it = []
    errorex = 1 + tolx
    erroref = 1 + tolf

    while True: # todo

        fx0 = 0 # todo
        d = 0 # todo
        
        x1 = 0 # todo
        fx1 = 0 # todo

        if x1 != 0:
            errorex = 0 # todo 
        else:
            errorex = 0 # todo
        
        erroref = 0 # todo
        
        x0 = x1
        it.append(x1)
        
    if len(it) == nmax:
        print('Corde : raggiunto massimo numero di iterazioni \n')
            
    return x1, it, len(it)
    
def newton(fname, fpname, x0, m = 1, nmax = 100, tolx = 1e-12, tolf = 1e-12):
        it = []
        errorex = 1 + tolx
        erroref = 1 + tolf

        while True: #todo
            fx0 = fname(x0)
            if False: # todo
                print(" derivata prima nulla in x0")
                return None, None,None
            
            d = 0 # todo 
            x1 = x0 - m * d

            if x1 != 0:
                errorex = 0 # todo
            else:
                errorex = 0 # todo 

            erroref = numpy.abs(fname(x1))
            it.append(x1)
            x0 = x1
          
        if len(it) == nmax:
            print('Newton: raggiunto massimo numero di iterazioni \n')

        return x1, it, len(it)
    
def secanti(fname, xm1, x0, nmax = 100, tolx = 1e-12, tolf = 1e-12):
        it = []
        errorex = 1 + tolx
        erroref = 1 + tolf

        while True: # todo
            fxm1 = 0 # todo
            fx0 = 0 # todo 
            d = 0 # todo 

            x1 = 0 # todo
            
            fx1 = fname(x1)
            if x1!=0:
                errorex= 0 # todo 
            else:
                errorex = 0 # todo
            
            it.append(x1);
            erroref = 0 # todo 
            xm1 = 0 # todo 
            x0 = 0# todo
       
        if len(it) == nmax:
            print('Secanti: raggiunto massimo numero di iterazioni \n')
        
        return x1, it, len(it)
    
def stima_ordine(xk,iterazioni):
    #Vedi dispensa allegata per la spiegazione
    k = iterazioni - 4
    p = numpy.log(abs(xk[k + 2] - xk[k + 3]) / abs(xk[k + 1] - xk[k + 2])) / numpy.log(abs(xk[k + 1] - xk[k + 2]) / abs(xk[k] - xk[k + 1]));
    
    return p

#Soluzione di sistemi di equazioni non lineari
def newton_raphson(initial_guess, F_numerical, J_Numerical, max_it = 100, tolX = 1e-12, tolF = 1e-12):
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    errore=[]
    
    while True: # todo
        jx = 0 # todo
        
        if False: # todo
            print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
            return None, None, None
        
        fx = 0# todo
        fx = fx.squeeze() 
        s = 0 # todo 
        Xnew = 0 # todo
        normaXnew = numpy.linalg.norm(Xnew, 1)

        if normaXnew != 0:
            erroreX = 0 # todo 
        else:
            erroreX = 0 # todo 
        
        errore.append(erroreX)
        fxnew = 0 # todo
        erroreF= numpy.linalg.norm(fxnew.squeeze(),1)
        X = Xnew
    
    return X,errore, len(errore)

def newton_raphson_corde(initial_guess, F_numerical, J_Numerical, max_it = 100, tolX = 1e-12, tolF = 1e-12):
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    errore=[]
    
    while True: # todo
        if len(errore): # todo
            jx = 0 # todo
        
            if False: # todo
                print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
                return None, None,None
        
        fx = 0 # todo
        fx = fx.squeeze() 
        s = 0 # todo 
        Xnew = 0 # todo
        normaXnew = numpy.linalg.norm(Xnew, 1)

        if normaXnew !=0:
            erroreX = 0 # todo 
        else:
            erroreX = 0 # todo 
        
        errore.append(erroreX)
        fxnew = 0 # todo
        erroreF = numpy.linalg.norm(fxnew.squeeze(), 1)
        X = Xnew
    
    return X, errore, len(errore)


def newton_raphson_sham(initial_guess, update, F_numerical, J_Numerical, max_it = 100,  tolX = 1e-12, tolF = 1e-12):
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    errore = []
    
    while True: # todo
        if len(errore): # toto
            jx = 0 # todo
            if False: # todo
                print("La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo")
                return None, None, None
        
        fx = 0 # todo
        fx = fx.squeeze() 
        s = 0 # todo 
        Xnew = 0 # todo
        normaXnew = numpy.linalg.norm(Xnew, 1)

        if normaXnew != 0:
            erroreX = 0 # todo 
        else:
            erroreX= 0 # todo
        
        errore.append(erroreX)
        fxnew = 0 # todo
        erroreF = numpy.linalg.norm(fxnew.squeeze(), 1)
        X = Xnew
    
    return X, errore, len(errore)

#Minimo di una funzion enon lineare
def newton_raphson_minimo(initial_guess, grad_func, Hessian_func, tolX, tolF, max_iterations):
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolX
    erroreX = 1 + tolF
    errore = []
    
    while True: # todo
        Hx = 0 # todo
        if False: # todo
            print("La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo")
            return None, None, None

        gfx = 0 # todo 
        gfx = gfx.squeeze() 
        s = 0 # todo 
        Xnew = 0 # todo 
        normaXnew = numpy.linalg.norm(Xnew, 1)

        if normaXnew != 0:
            erroreX = 0# todo
        else:
            erroreX = 0 # todo
            
        errore.append(erroreX)
        gfxnew = 0 # todo 
        erroreF = numpy.linalg.norm(gfxnew.squeeze(), 1)
        X = Xnew
    
    return X, errore, len(errore)

#Metodi Iterativi basati sullo splitting della matrice: jacobi, gauss-Seidel - Gauss_seidel SOR
def jacobi(A, b, x0, max_it = 100, toll = 1e-12):
    errore = 1000
    d = numpy.diag(A)
    n = A.shape[0]
    invM = numpy.diag(1/d)
    E = 0 # todo 
    F = 0 # todo 
    N = 0 # todo 
    T = 0 # todo 
    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = 0 # todo 
    print("raggio spettrale jacobi", raggiospettrale)

    it=0
    er = []
    while len(er) <= max_it and errore >= toll:
        x = 0 # todo 
        errore = 0# todo 
        er.append(errore)
        x0 = x.copy()

    return x, er, len(er)

def gauss_seidel(A, b, x0,max_it = 100, toll = 1e-12):
    errore = 1000
    d = 0 # todo 
    D = 0 # todo 
    E = 0 # todo 
    F = 0 # todo 
    M = 0 # todo 
    N = 0 # todo 
    T = 0 # todo 
    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = 0 # todo 
    print("raggio spettrale Gauss-Seidel ",raggiospettrale)

    er=[]
    while True: # todo 
        x = 0 # todo 
        errore = 0 # todo 
        er.append(errore)
        x0 = x.copy()

    return x, er, len(er)

def gauss_seidel_sor(A, b, x0, omega, max_it = 100, toll = 1e-12):
    errore = 1000
    d = 0 # todo 
    D = 0 # todo 
    E = 0 # todo 
    F = 0 # todo  
    Momega = D + omega * E
    Nomega = (1 - omega) * D - omega * F
    T = 0 # todo
    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = 0 # todo
    print("raggio spettrale Gauss-Seidel SOR ", raggiospettrale)
    
    M = 0 # todo 
    N = 0 # todo 

    xold = x0.copy()
    xnew = x0.copy()
    er=[]
    while len(er) <= max_it and errore >= toll:
        xtilde = 0 # todo 
        xnew = 0 # todo 
        errore = 0 # todo 
        er.append(errore)
        xold = xnew.copy()

    return xnew, er, len(er)

#Metodi di Discesa
def steepestdescent(A, b, x0, max_it = 100, toll = 1e-12):
    n, m = A.shape
    if n != m:
        print("Matrice non quadrata")
        return [],[]
    
    # inizializzare le variabili necessarie
    x = 0 
    r = 0 # todo 
    p = 0 # todo 
    it = 0
    nb = numpy.linalg.norm(b)
    errore = numpy.linalg.norm(r)/nb
    vec_sol = []
    vec_sol.append(x.copy())
    vet_r = []
    vet_r.append(errore)

    # utilizzare il metodo del gradiente per trovare la soluzione
    while True: # todo 
        it = it + 1 
        Ap = 0 # todo 
        alpha = 0 # todo
        x = 0 # todo         

        vec_sol.append(x.copy())
        r = 0 # todo 
        errore = numpy.linalg.norm(r) / nb
        vet_r.append(errore)
        p = 0 # todo 
        
    iterates_array = np.vstack([arr.T for arr in vec_sol])
    return x, vet_r, iterates_array, it

def conjugate_gradient(A, b, x0, max_it = 100, toll = 1e-12):
    n, m = A.shape
    if n != m:
        print("Matrice non quadrata")
        return [],[]
    
    # inizializzare le variabili necessarie
    x = x0
    
    r = 0 # todo 
    p = 0 # todo 
    it = 0
    nb = numpy.linalg.norm(b)
    errore = numpy.linalg.norm(r) / nb
    vec_sol=[]
    vec_sol.append(x0.copy())
    vet_r = []
    vet_r.append(errore)
    
    # utilizzare il metodo del gradiente coniugato per calcolare la soluzione
    while True: # todo 
        it=it+1
        Ap = 0 # todo A.dot(p)
        alpha = 0# todo 
        x = 0 # todo 
        vec_sol.append(x.copy())
        rtr_old = 0 # todo
        r = 0 # todo 
        gamma = 0 # todo 
        errore = numpy.linalg.norm(r)/nb
        vet_r.append(errore)
        p = 0 # todo 
   
    iterates_array = np.vstack([arr.T for arr in vec_sol])
    return x, vet_r, iterates_array, it

#Soluzione di sistemi sovradeterminati
def eqnorm(A,b):
    G = 0 # todo  
    f = 0 # todo 
    
    L = 0 # todo
    U = 0 # todo
    return x

def qrLS(A,b):
    n, m = A.shape[1]  # numero di colonne di A
    Q, R = scipy.linalg.qr(A)
    h = 0 # todo 
    x, _ = 0 # todo
    residuo = 0 # todo 
    return x, residuo

def SVDLS(A,b):
    m, n = A.shape  #numero di righe e  numero di colonne di A
    U, s, VT = scipy.linalg.svd(A)  
    
    V = VT.T
    thresh = numpy.spacing(1) * m * s[0] ##Calcolo del rango della matrice, numero dei valori singolari maggiori di una soglia
    k = 0 # todo 
    
    d = 0 # todo 
    d1 = 0 # todo 
    s1 = 0 # todo 
    
    c = 0 # todo 
    x = 0 # todo 
    residuo = 0 # todo 
    return x, residuo

#-----------Interpolazione
def plagr(xnodi,j):
    xzeri = numpy.zeros_like(xnodi)
    n = xnodi.size

    if j == 0:
       xzeri == 0 # todo 
    else:
       xzeri = numpy.append() # todo
    
    num = 0 
    den = 0 
    
    p = 0
    
    return p

def InterpL(x, y, xx):     
    n = 0 # todo
    m = 0 # todo 
    
    L = numpy.zeros((m,n))
    
    for j in 0: # todo
        p = 0 # todo 
        L[:, j] = 0 # todo 

    return 0 # todo 
