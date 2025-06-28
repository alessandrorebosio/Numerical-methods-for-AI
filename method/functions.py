# Zeri di funzione

import numpy.linalg
import math, numpy, scipy
import scipy.linalg


def sign(x):
    """
    Funzione segno che restituisce 1 se x è positivo, 0 se x è zero e -1 se x è negativo.
    """
    return math.copysign(1, x)


def metodo_bisezione(fname, a, b, tolX=1e-12):
    if sign(fname(a) * fname(b)) < 0:
        print("Non è possibile applicare il metodo di bisezione \n")
        return None, None, None

    it = []
    while math.abs(a - b) > tolX:
        c = (a + b) / 2
        it.append(c)

        if fname(c) == 0:
            return c, it, len(it)
        elif sign(fname(a) * fname(c)) < 0:
            b = c
        else:
            a = c

    return c, it, len(it)


def falsa_posizione(fname, a, b, max_it=100, tolX=1e-12):
    if sign(fname(a) * fname(b)) < 0:
        print("Metodo di bisezione non applicabile")
        return None, None, None

    erroreX = 1 + tolX
    prec = a
    it = []
    while len(it) < max_it and erroreX >= tolX:
        c = a - fname(a) * ((b - a) / (fname(a) - fname(b)))
        it.append(c)

        if fname(c) == 0:
            return c, it, len(it)
        elif sign(fname(a) * fname(c)) < 0:
            b = c
        else:
            a = c

        if c != 0:
            erroreX = math.abs(c - prec) / math.abs(c)  # todo
        else:
            erroreX = math.abs(c - prec)  # todo

        prec = c

    return c, it, len(it)


def corde(fname, coeff_ang, x0, max_it=100, tolX=1e-12, tolF=1e-12):
    # coeff_ang è il coefficiente angolare della retta che rimane fisso per tutte le iterazioni
    erroreX = 1 + tolX
    erroreF = 1 + tolF
    it = []
    while len(it) < max_it and erroreX >= tolX and erroreF >= tolF:
        d = fname(x0) / coeff_ang
        x1 = x0 - d

        if x1 != 0:
            erroreX = math.abs(d / x1)
        else:
            erroreX = math.abs(d)

        erroreF = math.abs(fname(x1))
        it.append(x1)
        x0 = x1

    if len(it) == max_it:
        print("Corde : raggiunto massimo numero di iterazioni \n")

    return x1, it, len(it)


def newton(fname, fpname, x0, m=1, max_it=100, tolX=1e-12, tolF=1e-12):
    erroreX = 1 + tolX
    erroreF = 1 + tolF
    it = []
    while len(it) < max_it and erroreX >= tolX and erroreF and tolF:
        if abs(fpname(x0) <= numpy.spacing(1)):
            print(" derivata prima nulla in x0")
            return None, None, None

        d = fname(x0) / fpname(x0)
        x1 = x0 - m * d

        if x1 != 0:
            erroreX = math.abs(d / x1)
        else:
            erroreX = math.abs(d)

        erroreF = numpy.abs(fname(x1))
        it.append(x1)
        x0 = x1

    if len(it) == max_it:
        print("Newton: raggiunto massimo numero di iterazioni \n")

    return x1, it, len(it)


def secanti(fname, xm1, x0, max_it=100, tolX=1e-12, tolF=1e-12):
    erroreX = 1 + tolX
    erroreF = 1 + tolF
    it = []
    while len(it) < max_it and erroreX >= tolX and erroreF >= tolF:
        d = fname(x0) * (x0 - xm1) / (fname(x0) - fname(xm1))
        x1 = x0 - d

        if x1 != 0:
            erroreX = math.abs(d / x1)
        else:
            erroreX = math.abs(d)

        erroreF = math.abs(fname(x1))
        it.append(x1)
        xm1 = x0
        x0 = x1

    if len(it) == max_it:
        print("Secanti: raggiunto massimo numero di iterazioni \n")

    return x1, it, len(it)


def stima_ordine(xk, iterazioni):
    k = iterazioni - 4
    return numpy.log(
        abs(xk[k + 2] - xk[k + 3]) / math.abs(xk[k + 1] - xk[k + 2])
    ) / numpy.log(math.abs(xk[k + 1] - xk[k + 2]) / math.abs(xk[k] - xk[k + 1]))


# Soluzione di sistemi di equazioni non lineari
def newton_raphson(
    initial_guess, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    it = []
    while len(it) < max_it and erroreX >= tolX and erroreF >= tolF:
        jx = J_Numerical(*X)
        if numpy.linalg.det(jx) == 0:
            print(
                "La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo"
            )
            return None, None, None

        fx = F_numerical(*X).squeeze()
        s = numpy.linalg.solve(jx, -fx)
        X = X + s
        norma = numpy.linalg.norm(X, 1)

        if norma != 0:
            erroreX = numpy.linalg.norm(s, 1) / norma
        else:
            erroreX = numpy.linalg.norm(s, 1)

        erroreF = numpy.linalg.norm(F_numerical(*X).squeeze(), 1)
        it.append(erroreX)

    return X, it, len(it)


def newton_raphson_corde(
    initial_guess, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    it = []
    while len(it) < max_it and erroreX >= tolX and erroreF >= tolF:
        if len(it) == 0:
            jx = J_Numerical(*X)
            if numpy.linalg.det(jx) == 0:
                print(
                    "La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo"
                )
                return None, None, None

        fx = F_numerical(*X).squeeze()
        s = numpy.linalg.solve(jx, -fx)
        X = X + s
        norma = numpy.linalg.norm(X, 1)

        if norma != 0:
            erroreX = numpy.linalg.norm(s, 1) / norma
        else:
            erroreX = numpy.linalg.norm(s, 1)

        erroreF = numpy.linalg.norm(F_numerical(*X).squeeze(), 1)
        it.append(erroreX)

    return X, it, len(it)


def newton_raphson_sham(
    initial_guess, update, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    it = []
    while len(it) < max_it and erroreX >= tolX and erroreF >= tolF:
        if len(it) % update == 0:
            jx = J_Numerical(*X)
            if numpy.linalg.det(jx) == 0:
                print(
                    "La matrice dello Jacobiano calcolata nell'iterato precedente non è a rango massimo"
                )
                return None, None, None

        fx = F_numerical(*X).squeeze()
        s = numpy.linalg.solve(jx, -fx)
        X = X + s
        norma = numpy.linalg.norm(X, 1)

        if norma != 0:
            erroreX = numpy.linalg.norm(s, 1) / norma
        else:
            erroreX = numpy.linalg.norm(s, 1)

        erroreF = numpy.linalg.norm(F_numerical(*X).squeeze(), 1)
        it.append(erroreX)

    return X, it, len(it)


# Minimo di una funzion enon lineare
def newton_raphson_minimo(
    initial_guess,
    grad_func,
    Hessian_func,
    max_it=100,
    tolX=1e-12,
    tolF=1e-12,
):
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolX
    erroreX = 1 + tolF
    it = []

    while len(it) < max_it and erroreX >= tolX and erroreF >= tolF:
        Hx = Hessian_func(*X)
        if numpy.linalg.det(Hx) == 0:
            print(
                "La matrice Hessiana calcolata nell'iterato precedente non è a rango massimo"
            )
            return None, None, None

        gfx = grad_func(*X).squeeze()
        s = numpy.linalg.solve(Hx, -gfx)
        X = X + s
        norma = numpy.linalg.norm(X, 1)

        if norma != 0:
            erroreX = numpy.linalg.norm(s, 1) / norma
        else:
            erroreX = numpy.linalg.norm(s, 1)

        erroreF = numpy.linalg.norm(grad_func(*X).squeeze(), 1)
        it.append(erroreX)

    return X, it, len(it)


# Metodi Iterativi basati sullo splitting della matrice: jacobi, gauss-Seidel - Gauss_seidel SOR
def jacobi(A, b, x0, max_it=100, toll=1e-12):
    errore = 1000
    d = numpy.diag(A)
    n = A.shape[0]
    invM = numpy.diag(1 / d)
    E = 0  # todo
    F = 0  # todo
    N = 0  # todo
    T = 0  # todo
    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = 0  # todo
    print("raggio spettrale jacobi", raggiospettrale)

    it = 0
    er = []
    while len(er) <= max_it and errore >= toll:
        x = 0  # todo
        errore = 0  # todo
        er.append(errore)
        x0 = x.copy()

    return x, er, len(er)


def gauss_seidel(A, b, x0, max_it=100, toll=1e-12):
    errore = 1000
    d = 0  # todo
    D = 0  # todo
    E = 0  # todo
    F = 0  # todo
    M = 0  # todo
    N = 0  # todo
    T = 0  # todo
    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = 0  # todo
    print("raggio spettrale Gauss-Seidel ", raggiospettrale)

    er = []
    while True:  # todo
        x = 0  # todo
        errore = 0  # todo
        er.append(errore)
        x0 = x.copy()

    return x, er, len(er)


def gauss_seidel_sor(A, b, x0, omega, max_it=100, toll=1e-12):
    errore = 1000
    d = 0  # todo
    D = 0  # todo
    E = 0  # todo
    F = 0  # todo
    Momega = D + omega * E
    Nomega = (1 - omega) * D - omega * F
    T = 0  # todo
    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = 0  # todo
    print("raggio spettrale Gauss-Seidel SOR ", raggiospettrale)

    M = 0  # todo
    N = 0  # todo

    xold = x0.copy()
    xnew = x0.copy()
    er = []
    while len(er) <= max_it and errore >= toll:
        xtilde = 0  # todo
        xnew = 0  # todo
        errore = 0  # todo
        er.append(errore)
        xold = xnew.copy()

    return xnew, er, len(er)


# Metodi di Discesa
def steepestdescent(A, b, x0, max_it=100, toll=1e-12):
    n, m = A.shape
    if n != m:
        print("Matrice non quadrata")
        return [], []

    # inizializzare le variabili necessarie
    x = 0
    r = 0  # todo
    p = 0  # todo
    it = 0
    nb = numpy.linalg.norm(b)
    errore = numpy.linalg.norm(r) / nb
    vec_sol = []
    vec_sol.append(x.copy())
    vet_r = []
    vet_r.append(errore)

    # utilizzare il metodo del gradiente per trovare la soluzione
    while True:  # todo
        it = it + 1
        Ap = 0  # todo
        alpha = 0  # todo
        x = 0  # todo

        vec_sol.append(x.copy())
        r = 0  # todo
        errore = numpy.linalg.norm(r) / nb
        vet_r.append(errore)
        p = 0  # todo

    iterates_array = np.vstack([arr.T for arr in vec_sol])
    return x, vet_r, iterates_array, it


def conjugate_gradient(A, b, x0, max_it=100, toll=1e-12):
    n, m = A.shape
    if n != m:
        print("Matrice non quadrata")
        return [], []

    # inizializzare le variabili necessarie
    x = x0

    r = 0  # todo
    p = 0  # todo
    it = 0
    nb = numpy.linalg.norm(b)
    errore = numpy.linalg.norm(r) / nb
    vec_sol = []
    vec_sol.append(x0.copy())
    vet_r = []
    vet_r.append(errore)

    # utilizzare il metodo del gradiente coniugato per calcolare la soluzione
    while True:  # todo
        it = it + 1
        Ap = 0  # todo A.dot(p)
        alpha = 0  # todo
        x = 0  # todo
        vec_sol.append(x.copy())
        rtr_old = 0  # todo
        r = 0  # todo
        gamma = 0  # todo
        errore = numpy.linalg.norm(r) / nb
        vet_r.append(errore)
        p = 0  # todo

    iterates_array = np.vstack([arr.T for arr in vec_sol])
    return x, vet_r, iterates_array, it


# Soluzione di sistemi sovradeterminati
def eqnorm(A, b):
    G = 0  # todo
    f = 0  # todo

    L = 0  # todo
    U = 0  # todo
    return x


def qrLS(A, b):
    n, m = A.shape[1]  # numero di colonne di A
    Q, R = scipy.linalg.qr(A)
    h = 0  # todo
    x, _ = 0  # todo
    residuo = 0  # todo
    return x, residuo


def SVDLS(A, b):
    m, n = A.shape  # numero di righe e  numero di colonne di A
    U, s, VT = scipy.linalg.svd(A)

    V = VT.T
    thresh = (
        numpy.spacing(1) * m * s[0]
    )  ##Calcolo del rango della matrice, numero dei valori singolari maggiori di una soglia
    k = 0  # todo

    d = 0  # todo
    d1 = 0  # todo
    s1 = 0  # todo

    c = 0  # todo
    x = 0  # todo
    residuo = 0  # todo
    return x, residuo


# -----------Interpolazione
def plagr(xnodi, j):
    xzeri = numpy.zeros_like(xnodi)
    n = xnodi.size

    if j == 0:
        xzeri == 0  # todo
    else:
        xzeri = numpy.append()  # todo

    num = 0
    den = 0

    p = 0

    return p


def InterpL(x, y, xx):
    n = 0  # todo
    m = 0  # todo

    L = numpy.zeros((m, n))

    for j in 0:  # todo
        p = 0  # todo
        L[:, j] = 0  # todo

    return 0  # todo
