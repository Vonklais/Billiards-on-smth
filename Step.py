import mpmath as mp

mp.dps = 50  # высокая точность

# Заданные корни (в порядке как на рисунке)
roots = [0, 15, 35, 40, 75]
z1, z2, z3, z4, z5 = roots
g = 2

# Полином и корень
def P(z):
    return (z-z1)*(z-z2)*(z-z3)*(z-z4)*(z-z5)



def sqrtP(z, branch=1):
    """Корень с выбором ветви (branch=1 для верхнего, -1 для нижнего листа)"""
    return branch * mp.sqrt(P(z))

def arc_integral(center, radius, func_upper, func_lower):
    upper = mp.quad(lambda t: func_upper(center + radius * mp.exp(1j * t)) * 1j * radius * mp.exp(1j * t),
                    [0, mp.pi])
    lower = mp.quad(lambda t: func_lower(center + radius * mp.exp(1j * t)) * 1j * radius * mp.exp(1j * t),
                    [mp.pi, 2*mp.pi])
    return upper + lower

print(mp.quad(lambda t: t+mp.j,
                    [2, 0]))