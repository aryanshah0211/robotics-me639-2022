from tarfile import SYMTYPE
import math
import sympy as sy
import numpy as np

n = 2 
#the number of joints for the 2R manipulator


q1 = sy.Symbol('q1')
q1_dot = sy.Symbol('q1_dot')
q1_doubledot = sy.Symbol('q1_doubledot')

q2 = sy.Symbol('q2')
q2_dot = sy.Symbol('q2_dot')
q2_doubledot = sy.Symbol('q2_doubledot')

m1 = sy.Symbol('m1')
m2 = sy.Symbol('m2')
l1 = sy.Symbol('l1')
l2 = sy.Symbol('l2')
g = 9.81

# things will change according to increase in link number 
D=np.array([[(m1*l1**2)/3+m2*l1**2, m2*l1*l2/2*sy.cos(q2-q1)],[m2*l1*l2/2*sy.cos(q2-q1), (m2*l2**2)/3]])
V=m1*g*l1/2*sy.sin(q1)+m2*g*l1*sy.sin(q1)+ m2*g*l2/2*sy.sin(q2)

#the potential energy derivative term
phi=np.array([[sy.diff(V, q1)],[sy.diff(V, q2)]])

q=np.array([[q1],[q2]])
q_dot=np.array([[q1_dot],[q2_dot]])
q_doubledot=np.array([[q1_doubledot],[q2_doubledot]])

#initialising the chirstoffel symbols of first kind
c=[0]*n

for k in range(n):
    for i in range(n):
        for j in range(n):
             c[k]+=0.5*(sy.diff(D[k][j], "q"+str(i+1))+sy.diff(D[k][i], "q"+str(j+1))-sy.diff(D[i][j], "q"+str(k+1)))*sy.Symbol("q"+str(i+1)+"_dot")*sy.Symbol("q"+str(j+1)+"_dot")

#final equation 
tau= D@q_doubledot + phi + np.transpose([c])

#the array will have the value of tau individual being tau1,tau2 for this case
print(tau)

