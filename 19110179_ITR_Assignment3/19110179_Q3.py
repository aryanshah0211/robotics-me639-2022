import math 
from math import pi
import numpy as np


#for the number of links
n=int(input("number of links: "))

#for the dh parameters
dh=[]
revolute_prismatic = [0]*n
for i in range (n):

    a=float(input("a_i: "))
    alpha=float(input("alpha_i (degree): "))*pi/180
    theta=float(input("theta_i (degree): "))*pi/180
    d=float(input("d_i: "))
    #now for prismatic or revolute
    revolute_prismatic[i] =  int(input("enter 0 if revolute and 1 if prismatic:"))
    dh.append([a,alpha,theta,d])

#the known dh matrix
A = []
for i in range (n):
    a,alpha,theta,d = dh[i]
    A.append((np.array([[np.cos(theta), -np.sin(theta)*np.cos(alpha), np.sin(theta)*np.sin(alpha), a*np.cos(theta)],[np.sin(theta), np.cos(theta)*np.cos(alpha), -np.cos(theta)*np.sin(alpha), a*np.sin(theta)],[0, np.sin(alpha), np.cos(alpha), d],[0, 0, 0, 1]])))

transformation = np.identity(4)

#for transformation matrix
for i in range (n):
    transformation = (transformation) @ (A[i])

print("the transformation matrix is:")
print(transformation)
print("position of end effector:")
print(transformation[0:3,3]) 
print("orientation of end effector:")
print(transformation[0:3,0:3])

#jacobian calculations
#initialising o and z matrices for jacobian calculation
o=[0]*(n+1)
o[0]=np.array([0, 0, 0])
z=[0]*(n+1)
z[0]=np.array([[0],[0],[1]])

for i in range (n):
    o[i+1]=transformation[0:3,3].T
    z[i+1]=transformation[0:3, 0:3]@z[0]

#calculates jacobian for every joint by downways concatenating (axis = 0) the Jv and Jw and adds it onto the jacobian_individual matrix
jacobian_individual = [0]*n
for i in range(n):
    if revolute_prismatic[i]==0:
        #so if revolute 
        jacobian_individual[i]=np.concatenate((np.cross(z[i].T,(o[len(n)]-o[i]).T).T,z[i]),axis=0)
    else:
        #prismatic
        jacobian_individual[i]=np.concatenate((z[i],np.array([[0,0,0]]).T),axis=0)

#jacobian matrix is made by sideways (axis =1) concatenating the elements of the jacobian_individual matrix
jacobian = jacobian_individual[0]
for i in range(1,n):
    jacobian=np.concatenate((jacobian,jacobian_individual[i]),axis=1)

print("jacobian:")
print(jacobian)

#for the velocity of end effector
q_dot = [0]*n
for i in range (n):
    q_dot[i] = float(print("enter the joints velocity q"+str(i+1)+"dot:"))

print("velocity of end effector:")
print(jacobian@q_dot[0:3])









