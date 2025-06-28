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
    """
    Metodo di bisezione per la ricerca dello zero di una funzione.
    Restituisce l'approssimazione dello zero, la lista degli iterati e il numero di iterazioni.
    """
    if sign(fname(a) * fname(b)) >= 0:
        print("Non è possibile applicare il metodo di bisezione")
        return None, None, None

    xk = []
    while abs(a - b) > tolX:
        c = a + (b - a) / 2
        xk.append(c)

        if fname(c) == 0:
            return c, xk, len(xk)
        elif sign(fname(a) * fname(c)) < 0:
            b = c
        else:
            a = c

    return c, xk, len(xk)


def falsa_posizione(fname, a, b, max_it=100, tolX=1e-12, tolF=1e-12):
    """
    Metodo della falsa posizione per la ricerca dello zero di una funzione.
    Restituisce l'approssimazione dello zero, la lista degli iterati e il numero di iterazioni.
    """
    if sign(fname(a) * fname(b)) >= 0:
        print("Metodo di bisezione non applicabile")
        return None, None, None

    erroreX = 1 + tolX
    erroreF = 1 + tolF
    prec = a
    xk = []
    while len(xk) < max_it and erroreF >= tolF and erroreX >= tolX:
        c = a - fname(a) * ((b - a) / (fname(b) - fname(a)))
        xk.append(c)

        if fname(c) == 0:
            return c, xk, len(xk)
        if sign(fname(a) * fname(c)) < 0:
            b = c
        else:
            a = c

        if c != 0:
            erroreX = abs(c - prec) / abs(c)
        else:
            erroreX = abs(c - prec)

        erroreF = abs(fname(c))
        prec = c

    return c, xk, len(xk)


def corde(fname, coeff_ang, x0, max_it=100, tolX=1e-12, tolF=1e-12):
    """
    Metodo delle corde per la ricerca dello zero di una funzione.
    Restituisce l'approssimazione dello zero, la lista degli iterati e il numero di iterazioni.
    """
    # coeff_ang è il coefficiente angolare della retta che rimane fisso per tutte le iterazioni
    erroreX = 1 + tolX
    erroreF = 1 + tolF
    xk = []
    while len(xk) < max_it and erroreX >= tolX and erroreF >= tolF:
        d = fname(x0) / coeff_ang
        x1 = x0 - d

        if x1 != 0:
            erroreX = abs(d / x1)
        else:
            erroreX = abs(d)

        erroreF = abs(fname(x1))
        xk.append(x1)
        x0 = x1

    if len(xk) == max_it:
        print("Corde : raggiunto massimo numero di iterazioni \n")

    return x1, xk, len(xk)


def newton(fname, fpname, x0, m=1, max_it=100, tolX=1e-12, tolF=1e-12):
    """
    Metodo di Newton per la ricerca dello zero di una funzione.
    Restituisce l'approssimazione dello zero, la lista degli iterati e il numero di iterazioni.
    """
    erroreX = 1 + tolX
    erroreF = 1 + tolF
    xk = []
    while len(xk) < max_it and erroreF >= tolF and erroreX >= tolX:
        if abs(fpname(x0)) <= numpy.spacing(1):
            print(" derivata prima nulla in x0")
            return None, None, None

        d = fname(x0) / fpname(x0)
        x1 = x0 - m * d

        if x1 != 0:
            erroreX = abs(d / x1)
        else:
            erroreX = abs(d)

        erroreF = numpy.abs(fname(x1))
        xk.append(x1)
        x0 = x1

    if len(xk) == max_it:
        print("Newton: raggiunto massimo numero di iterazioni \n")

    return x1, xk, len(xk)


def secanti(fname, xm1, x0, max_it=100, tolX=1e-12, tolF=1e-12):
    """
    Metodo delle secanti per la ricerca dello zero di una funzione.
    Restituisce l'approssimazione dello zero, la lista degli iterati e il numero di iterazioni.
    """
    erroreX = 1 + tolX
    erroreF = 1 + tolF
    xk = []
    while len(xk) < max_it and erroreX >= tolX and erroreF >= tolF:
        d = fname(x0) * (x0 - xm1) / (fname(x0) - fname(xm1))
        x1 = x0 - d

        if x1 != 0:
            erroreX = abs(d / x1)
        else:
            erroreX = abs(d)

        erroreF = abs(fname(x1))
        xk.append(x1)
        xm1 = x0
        x0 = x1

    if len(xk) == max_it:
        print("Secanti: raggiunto massimo numero di iterazioni \n")

    return x1, xk, len(xk)


def stima_ordine(xk, iterazioni):
    k = iterazioni - 4
    return numpy.log(
        abs(xk[k + 2] - xk[k + 3]) / abs(xk[k + 1] - xk[k + 2])
    ) / numpy.log(abs(xk[k + 1] - xk[k + 2]) / abs(xk[k] - xk[k + 1]))


# Soluzione di sistemi di equazioni non lineari
def newton_raphson(
    initial_guess, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    """
    Metodo di Newton-Raphson per la risoluzione di sistemi di equazioni non lineari.
    Restituisce la soluzione approssimata, la lista degli errori e il numero di iterazioni.
    """
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    er = []
    while len(er) < max_it and erroreX >= tolX and erroreF >= tolF:
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
        er.append(erroreX)

    return X, er, len(er)


def newton_raphson_corde(
    initial_guess, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    """
    Variante del metodo di Newton-Raphson che utilizza un'approccio a corde.
    Restituisce la soluzione approssimata, la lista degli errori e il numero di iterazioni.
    """
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    er = []
    while len(er) < max_it and erroreX >= tolX and erroreF >= tolF:
        if len(er) == 0:
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
        er.append(erroreX)

    return X, er, len(er)


def newton_raphson_sham(
    initial_guess, update, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    """
    Variante del metodo di Newton-Raphson che prevede un aggiornamento periodico dello Jacobiano.
    Restituisce la soluzione approssimata, la lista degli errori e il numero di iterazioni.
    """
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    er = []
    while len(er) < max_it and erroreX >= tolX and erroreF >= tolF:
        if len(er) % update == 0:
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
        er.append(erroreX)

    return X, er, len(er)


# Minimo di una funzion enon lineare
def newton_raphson_minimo(
    initial_guess,
    grad_func,
    Hessian_func,
    max_it=100,
    tolX=1e-12,
    tolF=1e-12,
):
    """
    Metodo di Newton-Raphson per la ricerca del minimo di una funzione.
    Restituisce la soluzione approssimata, la lista degli errori e il numero di iterazioni.
    """
    X = numpy.array(initial_guess, dtype=float)
    erroreF = 1 + tolX
    erroreX = 1 + tolF
    er = []
    while len(er) < max_it and erroreX >= tolX and erroreF >= tolF:
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
        er.append(erroreX)

    return X, er, len(er)


# Metodi Iterativi basati sullo splitting della matrice: jacobi, gauss-Seidel - Gauss_seidel SOR
def jacobi(A, b, x0, max_it=100, toll=1e-12):
    """
    Metodo iterativo di Jacobi per la risoluzione di sistemi lineari Ax = b.
    Restituisce la soluzione approssimata, la lista degli errori e il numero di iterazioni.
    """
    errore = 1000
    d = numpy.diag(A)
    n = A.shape[0]
    invM = numpy.diag(1 / d)

    N = -(numpy.tril(A, -1) + numpy.triu(A, 1))
    T = numpy.dot(invM, N)

    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = numpy.max(numpy.abs(autovalori))
    print("raggio spettrale jacobi", raggiospettrale)

    er = []
    while len(er) <= max_it and errore >= toll:
        x = (b + N @ x0) / d.reshape(n, 1)
        errore = numpy.linalg.norm(x - x0) / numpy.linalg.norm(x)
        er.append(errore)
        x0 = x.copy()

    return x, er, len(er)


def Lsolve(L, b):
    """
    Risoluzione con procedura forward di Lx=b con L triangolare inferiore.
    Restituisce la soluzione x e un flag di errore.
    """
    m, n = L.shape
    if n != m or numpy.all(numpy.diag(L)) != True:
        return [], 1

    x = numpy.zeros((n, 1))
    for i in range(n):
        s = numpy.dot(L[i, :i], x[:i])
        x[i] = (b[i] - s) / L[i, i]

    return x, 0


def Usolve(U, b):
    """
    Risoluzione con procedura backward di Ux=b con U triangolare superiore
     Input: U matrice triangolare superiore
            b termine noto
    Output: x: soluzione del sistema lineare
            flag=  0, se sono soddisfatti i test di applicabilità
                   1, se non sono soddisfatti

    """
    m, n = U.shape
    if n != m or numpy.all(numpy.diag(U)) != True:
        print("errore: matrice non quadrata")
        return [], 0

    x = numpy.zeros((n, 1))

    for i in range(n - 1, -1, -1):
        s = numpy.dot(U[i, i + 1 : n], x[i + 1 : n])
        x[i] = (b[i] - s) / U[i, i]

    return x, 1


def gauss_seidel(A, b, x0, max_it=100, toll=1e-12):
    """
    Metodo iterativo di Gauss-Seidel per la risoluzione di sistemi lineari Ax = b.
    Restituisce la soluzione approssimata, la lista degli errori e il numero di iterazioni.
    """
    errore = 1000
    d = numpy.diag(A)
    D = numpy.diag(d)
    E = numpy.tril(A, -1)
    F = numpy.triu(A, 1)
    M = D + E
    N = -F
    invM = numpy.linalg.inv(M)
    T = invM @ N
    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = numpy.max(numpy.abs(autovalori))
    print("raggio spettrale Gauss-Seidel ", raggiospettrale)

    er = []
    while len(er) <= max_it and errore >= toll:
        x, flag = Lsolve(M, b - F @ x0)
        errore = numpy.linalg.norm(x - x0) / numpy.linalg.norm(x)
        er.append(errore)
        x0 = x.copy()

    return x, er, len(er)


def gauss_seidel_sor(A, b, x0, omega, max_it=100, toll=1e-12):
    """
    Metodo iterativo di Gauss-Seidel SOR (Successive Over-Relaxation) per la risoluzione di sistemi lineari Ax = b.
    Permette di accelerare la convergenza tramite il parametro omega.
    Restituisce la soluzione approssimata, la lista degli errori e il numero di iterazioni.
    """
    errore = 1000
    d = numpy.diag(A)
    D = numpy.diag(d)
    E = numpy.tril(A, -1)
    F = numpy.triu(A, 1)
    Momega = D + omega * E
    Nomega = (1 - omega) * D - omega * F
    T = numpy.dot(numpy.linalg.inv(Momega), Nomega)
    autovalori = numpy.linalg.eigvals(T)
    raggiospettrale = numpy.max(numpy.abs(autovalori))
    M = D + E
    N = -F
    print("raggio spettrale Gauss-Seidel SOR ", raggiospettrale)

    xold = x0.copy()
    xnew = x0.copy()
    er = []
    while len(er) <= max_it and errore >= toll:
        xtilde, flag = Lsolve(M, b - numpy.dot(F, xold))
        xnew = (1 - omega) * xold + omega * xtilde
        errore = numpy.linalg.norm(xnew - xold) / numpy.linalg.norm(xnew)
        er.append(errore)
        xold = xnew.copy()

    return xnew, er, len(er)


# Metodi di Discesa
def steepestdescent(A, b, x0, max_it=100, toll=1e-12):
    """
    Metodo del gradiente (steepest descent) per la risoluzione di sistemi lineari Ax = b.
    Restituisce la soluzione approssimata, la lista degli errori, la sequenza degli iterati e il numero di iterazioni.
    """
    n, m = A.shape
    if n != m:
        print("Matrice non quadrata")
        return [], []

    x = x0.copy()
    r = b - A @ x
    it = 0
    nb = numpy.linalg.norm(b)
    errore = numpy.linalg.norm(r) / nb
    vec_sol = [x.copy()]
    vet_r = [errore]

    while errore >= toll and it < max_it:
        Ap = A @ r
        alpha = (r.T @ r) / (r.T @ Ap)
        x = x + alpha * r
        r = b - A @ x
        errore = numpy.linalg.norm(r) / nb
        vec_sol.append(x.copy())
        vet_r.append(errore)
        it += 1

    iterates_array = numpy.vstack([arr.T for arr in vec_sol])
    return x, vet_r, iterates_array, it


def conjugate_gradient(A, b, x0, max_it=100, toll=1e-12):
    """
    Metodo del gradiente coniugato per la risoluzione di sistemi lineari simmetrici e definiti positivi Ax = b.
    Restituisce la soluzione approssimata, la lista degli errori, la sequenza degli iterati e il numero di iterazioni.
    """
    n, m = A.shape
    if n != m:
        print("Matrice non quadrata")
        return [], []

    x = x0.copy()
    r = b - A @ x
    p = r.copy()
    it = 0
    nb = numpy.linalg.norm(b)
    errore = numpy.linalg.norm(r) / nb
    vec_sol = [x.copy()]
    vet_r = [errore]

    while errore >= toll and it < max_it:
        Ap = A @ p
        rtr = r.T @ r
        alpha = rtr / (p.T @ Ap)
        x = x + alpha * p
        r_new = r - alpha * Ap
        errore = numpy.linalg.norm(r_new) / nb
        vec_sol.append(x.copy())
        vet_r.append(errore)
        if errore < toll:
            break
        gamma = (r_new.T @ r_new) / rtr
        p = r_new + gamma * p
        r = r_new
        it += 1

    iterates_array = numpy.vstack([arr.T for arr in vec_sol])
    return x, vet_r, iterates_array, it


# Soluzione di sistemi sovradeterminati
def eqnorm(A, b):
    """
    Risoluzione di sistemi sovradeterminati Ax = b tramite il metodo dei minimi quadrati normali (normal equations).
    Restituisce la soluzione x.
    """
    G = A.T @ A
    f = A.T @ b

    L = scipy.linalg.cholesky(G, lower=True)
    U = L.T

    z, flag = Lsolve(L, f)
    if flag == 0:
        x, flag = Usolve(U, z)

    return x


def qrLS(A, b):
    """
    Risoluzione di sistemi sovradeterminati Ax = b tramite la decomposizione QR.
    Restituisce la soluzione x e il residuo.
    """
    n, m = A.shape
    Q, R = scipy.linalg.qr(A)
    h = Q.T @ b
    x, _ = Usolve(R[0:n, :], h[0:n])

    return x, numpy.linalg.norm(h[n:]) ** 2


def SVDLS(A, b):
    """
    Risoluzione di sistemi sovradeterminati Ax = b tramite la decomposizione ai valori singolari (SVD).
    Restituisce la soluzione x e il residuo.
    """
    m, n = A.shape
    U, s, VT = scipy.linalg.svd(A)

    V = VT.T
    # Calcolo del rango della matrice, numero dei valori singolari maggiori di una soglia
    thresh = numpy.spacing(1) * m * s[0]
    k = numpy.count_nonzero(s > thresh)

    d = U.T @ b
    d1 = d[:k].reshape(k, 1)
    s1 = s[:k].reshape(k, 1)

    c = d1 / s1
    x = V[:, :k] @ c

    return x, numpy.linalg.norm(d[k:]) ** 2


# -----------Interpolazione
def plagr(xnodi, j):
    """
    Calcola il polinomio fondamentale di Lagrange L_j(x) associato al nodo j.
    """
    xzeri = numpy.zeros_like(xnodi)
    n = xnodi.size

    if j == 0:
        xzeri = xnodi[1:n]
    else:
        xzeri = numpy.append(xnodi[0:j], xnodi[j + 1 : n])

    # Calcola i coefficienti del polinomio di grado n che si annulla nel vettore xzeri
    num = numpy.poly(xzeri)
    # Lo valuta nel nodo escluso (-jesimo)
    den = numpy.polyval(num, xnodi[j])

    return num / den


def InterpL(x, y, xx):
    """
    Interpolazione polinomiale: calcola il polinomio interpolante in xx dati i nodi x, y.
    """
    n = x.size
    m = xx.size

    L = numpy.zeros((m, n))
    for j in range(n):
        p = plagr(x, j)
        L[:, j] = numpy.polyval(p, xx)

    return L @ y
