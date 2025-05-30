import numpy as np


def A(a,b,c):
    return a, b, c,a,b,c,a,b

tre = [5,2]
weqr = [4]
q = A(*tre, *weqr)

w = q[:1]
e = q[1:]
print(w, e)