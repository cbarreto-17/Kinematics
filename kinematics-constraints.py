#!/usr/bin/env python3
import serial                     

#Combined kinematics to perform forward and inverse kinematics to find
#joint angles to desired location
import sys
import numpy as np
from datetime import datetime
from math import *
import random
## Forward Kinematics
def FK(angles,alpha,D,R):
    # Forward Kinematics
    # angles - joint angles (radians)
    # alpha - twist
    # d - offset
    # r - np.prod

    #Number of joints
    n=np.size(angles)-1
    p0=np.array([[0],[0],[0],[1]])

    # find np.cos and np.sin of angles and alpha
    st = np.array([])
    ct = np.array([])
    sa = np.array([])
    ca = np.array([])
    r = np.array([])
    d = np.array([])
    for i in range(0,n+1):
        st = np.concatenate((st,np.sin(angles[i])))
        ct = np.concatenate((ct,np.cos(angles[i])))
        sa = np.concatenate((sa,np.sin(alpha[i])))
        ca = np.concatenate((ca,np.cos(alpha[i])))
        r = np.concatenate((r,R[i]))
        d = np.concatenate((d,D[i]))


    # Find the homogenous transform from point 0 to unknown point
    H= (np.array([[ct[n],-st[n]*ca[n],st[n]*sa[n],r[n]*ct[n]],[st[n],ct[n]*ca[n],-ct[n]*sa[n],r[n]*st[n]],[0,sa[n],ca[n],d[n]],[0,0,0,1]]))  
    for i in range(n-1,-1,-1):
        Hi = (np.array([[ct[i],-st[i]*ca[i],st[i]*sa[i],r[i]*ct[i]],[st[i],ct[i]*ca[i],- ct[i]*sa[i],r[i]*st[i]],[0,sa[i],ca[i],d[i]],[0,0,0,1]]))
        H = np.matmul(Hi,H)

    point=np.matmul(H,p0)
    point=point[0:3]
    return point

## Jacobian
def Jacobian(angles,alpha,d,r):
    # Find the Jacobian matrix
    # angles - joint angles (radians)
    # alpha - twist
    # d - offset
    # r - np.prod

    #Number of joints
    n=np.prod(angles)

    # Assign variables to theta 
    t1 = angles[0]
    t2 = angles[1]
    t3 = angles[2]
    t4 = angles[3]
    t5 = angles[4]

    # Find the derivatives (NEED TO CHANGE THIS FOR VARYING DOFS ARM)
    j11 = -1*np.sin(t1)*(np.sin(t2+t3+t4)+cos(t2+t3)+np.cos(t2))
    j12 = -np.cos(t1)*(np.sin(t2+t3)-np.cos(t2+t3+t4)+np.sin(t2))
    j13 = np.cos(t1)*(np.cos(t2+t3+t4) - np.sin(t2 + t3))
    j14 = np.cos(t2 + t3 + t4)*np.cos(t1)
    j15 = [0]
    j21 = np.cos(t1)*(np.sin(t2+t3+t4)+np.cos(t2+t3)+np.cos(t2))
    j22 = -np.sin(t1)*(np.sin(t2+t3) - np.cos(t2+t3+t4) + np.sin(t2))
    j23 = np.sin(t1)*(np.cos(t2+t3+t4)-np.sin(t2+t3))
    j24 = np.cos(t2+t3+t4)*np.sin(t1)
    j25 = [0]
    j31 = [0]
    j32 = -np.sin(t2 + t3 + t4)-np.cos(t2+t3)-np.cos(t2)
    j33 = -np.sin(t2 + t3 + t4)-np.cos(t2+t3)
    j34 = -np.sin(t2+t3+t4)
    j35 = [0]

    #Finally, make the jacobian matrix
    r1 = np.concatenate((j11,j12,j13,j14,j15))
    r2 = np.concatenate((j21,j22,j23,j24,j25))
    r3 = np.concatenate((j31,j32,j33,j34,j35))
    J = np.concatenate((r1,r2,r3),axis = 0)
    J = np.reshape(J,(3,5))
    return J
    
## Find distance to a point
def Dist(point1,point2):

    x1=point1[0]
    y1=point1[1]
    z1=point1[2]
    x2=point2[0]
    y2=point2[1]
    z2=point2[2]
    d=sqrt(np.square(x2 - x1) + np.square(y2 - y1) + np.square(z2 - z1))
    return d

## Joint Limit constraints 
def jointlimits(angles):
    

    t1,t2,t3,t4,t5 = angles[0],angles[1],angles[2],angles[3],angles[4]


    # if t1 < 0 or t1 > pi:
    #     angles[0]  = np.deg2rad(np.floor(np.random.rand()*(180+1)))
    
    # if t2 < -5*pi/6 or t2 > pi/6:
    #     angles[1]  = np.deg2rad(np.floor(np.random.rand()*(180+1))-150)
    
    # if t3 < pi/3 or t3 > 4*pi/3:
    #     angles[1]  = np.deg2rad(np.floor(np.random.rand()*(180+1))+60)
    
    # if t4 < pi/2 or t4 > 3*pi/2:
    #     angles[1]  = np.deg2rad(np.floor(np.random.rand()*(180+1))+90)
    # if t5 < 0 or t5 > pi:
    #     angles[1]  = np.deg2rad(np.floor(np.random.rand()*(180+1)))
    if t1 < 0 or t1 > pi:
        angles[0]  = random.uniform(0, pi)
    
    if t2 < -5*pi/6 or t2 > pi/6:
        angles[1]  = random.uniform(-5*pi/6, pi/6)
    
    if t3 < pi/3 or t3 > 4*pi/3:
        angles[2]  = random.uniform(pi/3, 4*pi/3)
    
    if t4 < pi/2 or t4 > 4*pi/3:
        angles[3]  = random.uniform(pi/2, 4*pi/3)
    
    return angles


## Inverse Kinematics
def IK(angles,e,g,b,alpha,d,r):
    # Inverse Kinematics
    # angles - joint angles (radians)
    # alpha - twist
    # d - offset
    # r - np.prod
    # 0 < beta < 1.
    # g = desired position
    # e = current position

    #1. Compute J(q).
    J=Jacobian(angles,alpha,d,r)

    #2. Select an increment that will move e closer to g, where 0 < ? ? 1.
    de= b*(g - e)

    #3. Compute the change in the joints variables that will achieve the end-effector increment selected in the previous step
    JT = np.transpose(J)
    dangles=np.matmul(JT,de)

    #4. Apply this change to the joint variable vector, as q = q + ?q.
    angles=angles+dangles
    #print('Test angles (before)')
    #print(angles)
    angles = jointlimits(angles)

   # print('Test angles (after)')
    #print(angles)    

    #5. Update e by computing the forward kinematics of the manipulator with the updated q.
    point=FK(angles,alpha,d,r)
   # print('point')
   # print(point)
    return angles,point


def kinematics(xg,yg,zg):

    # No scientific notations
    np.set_printoptions(suppress=True)
    
    # Given desired location (x,y,g)
    giv_point=np.array([[xg],[yg],[zg]])
    pi = np.pi

    # Link sizes - L = [l1,l2,l3,l4]
    #l1 = 0.066675
    #l2 = 0.104775
    #l3 = 0.0889
    #l4 = 0.1778
    l1,l2,l3,l4=1,1,1,1
    L=np.array([l1,l2,l3,l4])

    #current joint angles (radians)
    angles=np.array([[0],[0],[0],[0],[0]])
    alpha=np.array([[-pi/2],[0],[0],[pi/2],[0]])
    r=np.array([[0],[L[1]],[L[2]],[0],[0]])
    d=np.array([[L[0]],[0],[0],[0],[L[3]]])


    #Decrease the tolerance to increase accuracy
    tolerance=0.001

    # Find current position based on the current joint angles
    point = FK(angles,alpha,d,r)
    error=np.inf
    b=0.1

    # kinematics.m:32
    # Iterate until joint angles found to the desired location
    while error > tolerance:

        angles,point=IK(angles,point,giv_point,b,alpha,d,r)
        error=Dist(point,giv_point)
        #print('Calculated position')
        #print(point)
    
    print('Desired Location')
    print(giv_point)
    print('Calculated position')
    print(point)
    print('Accuracy')
    print(tolerance)
    print('Angles (before)')
    print(np.rad2deg(angles))
    Z=np.array([[0],[5*pi/6],[-pi/3],[-pi/2],[0]])
    M = np.array([0,0,-1,0,0])
    angles= np.add(angles, Z)
    angles[2] = -1*angles[2]
    print('Angles (after)')
    print(np.rad2deg(angles))

    angles_actual = str().join(','.join('%0.3f' %y for y in np.rad2deg(angles)))
    ser = serial.Serial('/dev/ttyUSB0', 9600)
    ser.write(bytearray(angles_actual, 'utf8'))
    print(bytearray(angles_actual, 'utf8'))
    ser.close()
# Time taken
startTime = datetime.now()

# Call the kinematics function   
#kinematics(-0.17002123056796467,0.2241048910901859,0.22375420710815683)      #Claw points straight up (0,0,4)
#kinematics(0,2,2)  
print('TIME TAKEN:')
print(datetime.now() - startTime)
