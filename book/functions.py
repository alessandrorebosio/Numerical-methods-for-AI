import math, numpy as np, scipy
import scipy.linalg

# Ax = b con m = n

# piccole dimensioni e densa:
# ben condizionata: Gauss
# mal condizionata: QR
# simmetrica definita positiva: Choleski

# grandi dimensioni e sparsa:
# diagonale strettamente dominante: jacobi, Gauss
# simmetrica def. positiva: Gauss o metodi di discesa

# Ax = b con m > n

# ben condizionata e rango massimo: Equazioni normali
# mal condizionata e rango massimo: QRLS
# mal condizionata e non ha rango massimo: SVDLS


def Usolve(U, b):
    """
    Risoluzione con procedura backward di Ux=b con U triangolare superiore
     Input: U matrice triangolare superiore
            b termine noto
    Output: x: soluzione del sistema lineare
            flag=  0, se sono soddisfatti i test di applicabilità
                   1, se non sono soddisfatti

    """
    # test dimensione
    m, n = U.shape
    flag = 0
    if n != m:
        print("errore: matrice non quadrata")
        flag = 1
        x = []
        return x, flag

    # Test singolarita'
    if np.all(np.diag(U)) != True:
        print("el. diag. nullo - matrice triangolare superiore")
        x = []
        flag = 1
        return x, flag
    # Preallocazione vettore soluzione
    x = np.zeros((n, 1))

    for i in range(n - 1, -1, -1):
        s = np.dot(
            U[i, i + 1 : n], x[i + 1 : n]
        )  # scalare=vettore riga * vettore colonna
        x[i] = (b[i] - s) / U[i, i]

    return x, flag


def Lsolve(L, b):
    """
    Risoluzione con procedura forward di Lx=b con L triangolare inferiore
     Input: L matrice triangolare inferiore
            b termine noto
    Output: x: soluzione del sistema lineare
            flag=  0, se sono soddisfatti i test di applicabilità
                   1, se non sono soddisfatti
    """
    # test dimensione
    m, n = L.shape
    flag = 0
    if n != m:
        print("errore: matrice non quadrata")
        flag = 1
        x = []
        return x, flag

    # Test singolarita'
    if np.all(np.diag(L)) != True:
        print("el. diag. nullo - matrice triangolare inferiore")
        x = []
        flag = 1
        return x, flag
    # Preallocazione vettore soluzione
    x = np.zeros((n, 1))

    for i in range(n):
        s = np.dot(L[i, :i], x[:i])  # scalare=vettore riga * vettore colonna
        x[i] = (b[i] - s) / L[i, i]

    return x, flag


def sign(x):
    """
    Funzione segno che restituisce 1 se x è positivo, 0 se x è zero e -1 se x è negativo.
    """
    return math.copysign(1, x)


# Zeri di funzione
def metodo_bisezione(fname, a, b, tolX=1e-12):
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
    if sign(fname(a) * fname(b)) >= 0:
        print("Falsa posizione non applicabile")
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


# coeff_ang = (f(b) - f(a)) / (b - a)
def corde(fname, coeff_ang, x0, max_it=100, tolX=1e-12, tolF=1e-12):
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
    erroreX = 1 + tolX
    erroreF = 1 + tolF
    xk = []
    while len(xk) < max_it and erroreF >= tolF and erroreX >= tolX:
        if abs(fpname(x0)) <= np.spacing(1):
            print("derivata prima nulla in x0")
            return None, None, None

        d = fname(x0) / fpname(x0)
        x1 = x0 - m * d

        if x1 != 0:
            erroreX = abs(d / x1)
        else:
            erroreX = abs(d)

        erroreF = np.abs(fname(x1))
        xk.append(x1)
        x0 = x1

    if len(xk) == max_it:
        print("Newton: raggiunto massimo numero di iterazioni \n")

    return x1, xk, len(xk)


def secanti(fname, xm1, x0, max_it=100, tolX=1e-12, tolF=1e-12):
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
    return np.log(abs(xk[k + 2] - xk[k + 3]) / abs(xk[k + 1] - xk[k + 2])) / np.log(
        abs(xk[k + 1] - xk[k + 2]) / abs(xk[k] - xk[k + 1])
    )


# Soluzione di sistemi di equazioni non lineari
def newton_raphson(
    initial_guess, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    X = np.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    er = []
    while len(er) < max_it and erroreX >= tolX and erroreF >= tolF:
        jx = J_Numerical(*X)
        if np.linalg.det(jx) == 0:
            print("La matrice dello Jacobiano non è a rango massimo")
            return None, None, None

        fx = F_numerical(*X).squeeze()
        s = np.linalg.solve(jx, -fx)
        X = X + s

        norma = np.linalg.norm(X, 1)
        if norma != 0:
            erroreX = np.linalg.norm(s, 1) / norma
        else:
            erroreX = np.linalg.norm(s, 1)

        erroreF = np.linalg.norm(F_numerical(*X).squeeze(), 1)
        er.append(erroreX)

    return X, er, len(er)


def newton_raphson_corde(
    initial_guess, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    X = np.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    er = []
    while len(er) < max_it and erroreX >= tolX and erroreF >= tolF:
        if len(er) == 0:
            jx = J_Numerical(*X)
            if np.linalg.det(jx) == 0:
                print("La matrice dello Jacobiano non è a rango massimo")
                return None, None, None

        fx = F_numerical(*X).squeeze()
        s = np.linalg.solve(jx, -fx)
        X = X + s

        norma = np.linalg.norm(X, 1)
        if norma != 0:
            erroreX = np.linalg.norm(s, 1) / norma
        else:
            erroreX = np.linalg.norm(s, 1)

        erroreF = np.linalg.norm(F_numerical(*X).squeeze(), 1)
        er.append(erroreX)

    return X, er, len(er)


def newton_raphson_sham(
    initial_guess, update, F_numerical, J_Numerical, max_it=100, tolX=1e-12, tolF=1e-12
):
    X = np.array(initial_guess, dtype=float)
    erroreF = 1 + tolF
    erroreX = 1 + tolX
    er = []
    while len(er) < max_it and erroreX >= tolX and erroreF >= tolF:
        if len(er) % update == 0:
            jx = J_Numerical(*X)
            if np.linalg.det(jx) == 0:
                print("La matrice dello Jacobiano non è a rango massimo")
                return None, None, None

        fx = F_numerical(*X).squeeze()
        s = np.linalg.solve(jx, -fx)
        X = X + s

        norma = np.linalg.norm(X, 1)
        if norma != 0:
            erroreX = np.linalg.norm(s, 1) / norma
        else:
            erroreX = np.linalg.norm(s, 1)

        erroreF = np.linalg.norm(F_numerical(*X).squeeze(), 1)
        er.append(erroreX)

    return X, er, len(er)


# Minimo di una funzione non lineare
def newton_raphson_minimo(
    initial_guess,
    grad_func,
    Hessian_func,
    max_it=100,
    tolX=1e-12,
    tolF=1e-12,
):
    X = np.array(initial_guess, dtype=float)
    erroreF = 1 + tolX
    erroreX = 1 + tolF
    er = []
    while len(er) < max_it and erroreX >= tolX and erroreF >= tolF:
        Hx = Hessian_func(*X)
        if np.linalg.det(Hx) == 0:
            print("La matrice Hessiana non è a rango massimo")
            return None, None, None

        gfx = grad_func(*X).squeeze()
        s = np.linalg.solve(Hx, -gfx)
        X = X + s
        norma = np.linalg.norm(X, 1)

        if norma != 0:
            erroreX = np.linalg.norm(s, 1) / norma
        else:
            erroreX = np.linalg.norm(s, 1)

        erroreF = np.linalg.norm(grad_func(*X).squeeze(), 1)
        er.append(erroreX)

    return X, er, len(er)


# Metodi Iterativi basati sullo splitting della matrice: jacobi, gauss-Seidel - Gauss_seidel SOR
def jacobi(A, b, x0, max_it=100, toll=1e-12):
    errore = 1000
    d = np.diag(A)
    n = A.shape[0]
    invM = np.diag(1 / d)

    N = -(np.tril(A, -1) + np.triu(A, 1))
    T = np.dot(invM, N)

    autovalori = np.linalg.eigvals(T)
    raggiospettrale = np.max(np.abs(autovalori))
    print("raggio spettrale jacobi", raggiospettrale)

    er = []
    while len(er) <= max_it and errore >= toll:
        x = (b + N @ x0) / d.reshape(n, 1)
        errore = np.linalg.norm(x - x0) / np.linalg.norm(x)
        er.append(errore)
        x0 = x.copy()

    return x, er, len(er)


# x0 = np.zeros((b.size, 1))
def gauss_seidel(A, b, x0, max_it=100, toll=1e-12):
    errore = 1000
    d = np.diag(A)
    D = np.diag(d)
    E = np.tril(A, -1)
    F = np.triu(A, 1)
    M = D + E
    N = -F
    invM = np.linalg.inv(M)
    T = invM @ N
    autovalori = np.linalg.eigvals(T)
    raggiospettrale = np.max(np.abs(autovalori))
    print("raggio spettrale Gauss-Seidel ", raggiospettrale)

    er = []
    while len(er) <= max_it and errore >= toll:
        x, flag = Lsolve(M, b - F @ x0)
        errore = np.linalg.norm(x - x0) / np.linalg.norm(x)
        er.append(errore)
        x0 = x.copy()

    return x, er, len(er)


def gauss_seidel_sor(A, b, x0, omega, max_it=100, toll=1e-12):
    errore = 1000
    d = np.diag(A)
    D = np.diag(d)
    E = np.tril(A, -1)
    F = np.triu(A, 1)
    Momega = D + omega * E
    Nomega = (1 - omega) * D - omega * F
    T = np.dot(np.linalg.inv(Momega), Nomega)
    autovalori = np.linalg.eigvals(T)
    raggiospettrale = np.max(np.abs(autovalori))
    M = D + E
    N = -F
    print("raggio spettrale Gauss-Seidel SOR ", raggiospettrale)

    xold = x0.copy()
    xnew = x0.copy()
    er = []
    while len(er) <= max_it and errore >= toll:
        xtilde, flag = Lsolve(M, b - np.dot(F, xold))
        xnew = (1 - omega) * xold + omega * xtilde
        errore = np.linalg.norm(xnew - xold) / np.linalg.norm(xnew)
        er.append(errore)
        xold = xnew.copy()

    return xnew, er, len(er)


# Metodi di Discesa
def steepestdescent(A, b, x0, max_it=100, toll=1e-12):
    n, m = A.shape
    if n != m:
        print("Matrice non quadrata")
        return [], []

    x = x0.copy()
    r = b - A @ x
    nb = np.linalg.norm(b)
    errore = np.linalg.norm(r) / nb
    vec_sol = [x.copy()]
    vet_r = [errore]
    it = 0

    while errore >= toll and it < max_it:
        Ap = A @ r
        alpha = (r.T @ r) / (r.T @ Ap)
        x = x + alpha * r
        r = b - A @ x
        errore = np.linalg.norm(r) / nb
        vec_sol.append(x.copy())
        vet_r.append(errore)
        it += 1

    iterates_array = np.vstack([arr.T for arr in vec_sol])
    return x, vet_r, iterates_array, it


def conjugate_gradient(A, b, x0, max_it=100, toll=1e-12):
    n, m = A.shape
    if n != m:
        print("Matrice non quadrata")
        return [], []

    x = x0.copy()
    r = b - A @ x
    p = r.copy()
    it = 0
    nb = np.linalg.norm(b)
    errore = np.linalg.norm(r) / nb
    vec_sol = [x.copy()]
    vet_r = [errore]

    while errore >= toll and it < max_it:
        it += 1
        Ap = A @ p
        rtr = r.T @ r
        alpha = rtr / (p.T @ Ap)
        x = x + alpha * p
        r_new = r - alpha * Ap
        errore = np.linalg.norm(r_new) / nb
        vec_sol.append(x.copy())
        vet_r.append(errore)
        if errore < toll:
            break
        gamma = (r_new.T @ r_new) / rtr
        p = r_new + gamma * p
        r = r_new

    iterates_array = np.vstack([arr.T for arr in vec_sol])
    return x, vet_r, iterates_array, it


# Soluzione di sistemi sovradeterminati
def eqnorm(A, b):
    G = A.T @ A
    condG = np.linalg.cond(G)
    print("Indice di condizionamento di G ", condG)
    f = A.T @ b

    L = scipy.linalg.cholesky(G, lower=True)
    U = L.T

    z, flag = Lsolve(L, f)
    if flag == 0:
        x, flag = Usolve(U, z)

    return x


def qrLS(A, b):
    n = A.shape[1]
    Q, R = scipy.linalg.qr(A)
    h = Q.T @ b
    x, _ = Usolve(R[0:n, :], h[0:n])

    return x, np.linalg.norm(h[n:]) ** 2


def SVDLS(A, b):
    m, n = A.shape
    U, s, VT = scipy.linalg.svd(A)

    V = VT.T

    thresh = np.spacing(1) * m * s[0]
    k = np.count_nonzero(s > thresh)
    print("rango =", k)

    d = U.T @ b
    d1 = d[:k].reshape(k, 1)
    s1 = s[:k].reshape(k, 1)

    c = d1 / s1
    x = V[:, :k] @ c

    return x, np.linalg.norm(d[k:]) ** 2


# Interpolazione
def plagr(xnodi, j):
    xzeri = np.zeros_like(xnodi)
    n = xnodi.size

    if j == 0:
        xzeri = xnodi[1:n]
    else:
        xzeri = np.append(xnodi[0:j], xnodi[j + 1 : n])

    num = np.poly(xzeri)
    den = np.polyval(num, xnodi[j])

    return num / den


def InterpL(x, y, xx):
    n = x.size
    m = xx.size

    L = np.zeros((m, n))
    for j in range(n):
        p = plagr(x, j)
        L[:, j] = np.polyval(p, xx)

    return L @ y


# Fattorizzazione di Gauss
def metodo_gauss(A, b):
    PT, L, U = scipy.linalg.lu(A)
    P = PT.T
    y, flag = Lsolve(L, P @ b)
    if not flag:
        x, flag = Usolve(U, y)

    return x


# Fattorizzazione QR
def metodo_householder(A, b):
    Q, R = scipy.linalg.qr(A)
    y = Q.T @ b
    x, _ = Usolve(R, y)

    return x


# Fattorizzazione di Choleski
def metodo_choleski(A, b):
    L = scipy.linalg.cholesky(A, lower=True)
    y, flag = Lsolve(L, b)
    if not flag:
        x, flag = Usolve(L.T, y)

    return x
