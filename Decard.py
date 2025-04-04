import math
import numpy as np

def Z_calc_3D(X, Y, a, b, c, sign = 1.0):
    return sign*c*math.sqrt(abs(1-(X/a)**2 - (Y/b)**2))

def Zx_3D(X, Y, Z, a, b, c):
    return -X*c*c/(Z*a*a)

def Zy_3D(X, Y, Z, a, b, c):
    return -Y*c*c/(Z*b*b)

def D_Z_calc_3D(X, Y, D_X, D_Y, a, b, c, sign = 1.0):
    Z = Z_calc_3D(X, Y, a, b, c, sign)
    return Zx_3D(X, Y, Z, a, b, c)*D_X + Zy_3D(X, Y, Z, a, b, c)*D_Y

def Decard_step_3D(Xold, Yold, D_Xold, D_Yold, a, b, c, signZ = 1.0):
    Zold = Z_calc_3D(Xold, Yold, a, b, c, signZ)
    D_Zold = D_Z_calc_3D(Xold, Yold, D_Xold, D_Yold, a, b, c, signZ)

    Xnew = Xold + D_Xold
    Ynew = Yold + D_Yold
    Znew = Z_calc_3D(Xnew, Ynew, a, b, c, signZ)#проблема со знаком
    if D_Zold>(Znew - Zold)/2:
        Znew = -Znew
        signZ = -signZ

    q = np.array([[Znew**2 * Zold**2 + (Xnew**2 * c**4 * Zold**2)/a**4, (Xnew * Ynew * c**4 * Zold**2)/(a**2 * b**2)],
                  [(Xnew * Ynew * c**4 * Zold**2)/(a**2 * b**2), Znew**2 * Zold**2 + (Ynew**2 * c**4 * Zold**2)/b**4]])
    w = np.array([c**4 * Znew**2 * (Xold*D_Xold/a**2 + Yold*D_Yold/b**2)*Xnew + D_Xold * Zold**2 * Znew**2,
                  c**4 * Znew**2 * (Xold*D_Xold/a**2 + Yold*D_Yold/b**2)*Ynew + D_Yold * Zold**2 * Znew**2])

    #ПРОВЕРИТЬ!!!
    D_Xnew, D_Ynew = np.linalg.solve(q,w)

    if D_Xnew == 0:
        D_Xnew = 1e-5
    if D_Ynew == 0:
        D_Ynew = 1e-5
    return Xnew, Ynew, D_Xnew, D_Ynew, signZ