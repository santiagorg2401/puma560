#!/usr/bin/env python3

import sympy
import rospy
import numpy as np 
import math as m
import sys

from sympy import cos, sin, pi
from rospy.exceptions import ROSInterruptException
from std_msgs.msg import Float64
from control_msgs.msg import JointControllerState

class puma560_control:
    def __init__(self):

        #Init ros node.
        rospy.init_node("puma560_control", anonymous=True)
        
        # Set up subscribers.
        # self.th1_sub = rospy.Subscriber("/puma560_description/j1_pc/state", JointControllerState, queue_size=1)
        # self.th2_sub = rospy.Subscriber("/puma560_description/j2_pc/state", JointControllerState, queue_size=1)
        # self.th3_sub = rospy.Subscriber("/puma560_description/j3_pc/state", JointControllerState, queue_size=1)
        # self.th4_sub = rospy.Subscriber("/puma560_description/j4_pc/state", JointControllerState, queue_size=1)
        # self.th5_sub = rospy.Subscriber("/puma560_description/j5_pc/state", JointControllerState, queue_size=1)
        # self.th6_sub = rospy.Subscriber("/puma560_description/j6_pc/state/process_value", JointControllerState, queue_size=1)

        # Set up publishers.
        self.th1_pub = rospy.Publisher("/puma560_description/j1_pc/command", Float64, queue_size=1)
        self.th2_pub = rospy.Publisher("/puma560_description/j2_pc/command", Float64, queue_size=1)
        self.th3_pub = rospy.Publisher("/puma560_description/j3_pc/command", Float64, queue_size=1)
        self.th4_pub = rospy.Publisher("/puma560_description/j4_pc/command", Float64, queue_size=1)
        self.th5_pub = rospy.Publisher("/puma560_description/j5_pc/command", Float64, queue_size=1)
        #self.th6_pub = rospy.Publisher("/puma560_description/j6_pc/command", Float64, queue_size=1)

        self.publisher(0,0,0,0,0)

        #Puma 560 DH parameters.
        a2, a3, d3, d4 = 0.4318, 0.0203, 0.1500, 0.4318
        self.a = [0, 0, a2, a3, 0, 0]
        self.d = [0, 0, d3, d4, 0, 0]

        # self.Toriginefector, self.b, self.T36 = self.sym_fwd_kinematics()

    def sym_fwd_kinematics(self):
        i = range(6)        #Number of links.
        T = []
        Toriginefector = 1

        thi, alj, aj, di = sympy.symbols("thi, alj, aj, di")

        Tij = sympy.Matrix([[cos(thi),          -sin(thi),          0,        aj], 
                            [sin(thi)*cos(alj), cos(thi)*cos(alj), -sin(alj), -sin(alj)*di], 
                            [sin(thi)*sin(alj), cos(thi)*sin(alj), cos(alj),  cos(alj)*di],
                            [0,                 0,                 0,         1      ]])

        #Puma 560 DH parameters.
        th1, th2, th3, th4, th5, th6, a2, a3, d3, d4 = sympy.symbols("th1, th2, th3, th4, th5, th6, a2, a3, d3, d4")
        al = [0, -pi/2, 0, -pi/2, pi/2, -pi/2]
        a = [0, 0, a2, a3, 0, 0]
        d = [0, 0, d3, d4, 0, 0]
        th = [th1, th2, th3, th4, th5, th6]

        for j in i:
            Tx = Tij.subs([(thi, th[j]), (alj, al[j]), (aj, a[j]), (di, d[j])])
            T.append(Tx)

        #Direct kinematics.
        for j in i:
            Toriginefector = Toriginefector*T[j]

        Toriginefector = sympy.simplify(Toriginefector)

        r11, r12, r13, r21, r22, r23, r31, r32, r33, px, py, pz = sympy.symbols("r11, r12, r13, r21, r22, r23, r31, r32, r33, px, py, pz")
        Tsym = sympy.Matrix([[r11, r12, r13, px],
                            [r21, r22, r23, py],
                            [r31, r32, r33, pz],
                            [0, 0, 0, 1]])

        T03inv = sympy.simplify((T[0]*T[1]*T[2]).inv())
        b = T03inv*Tsym
        T36 = sympy.simplify(sympy.simplify(T03inv*Toriginefector))

        return Toriginefector, b, T36

    def publisher(self, th1, th2, th3, th4, th5):
        th1_m, th2_m, th3_m, th4_m, th5_m = Float64(th1), Float64(th2), Float64(th3), Float64(th4), Float64(th5)
        self.th1_pub.publish(th1_m)
        self.th2_pub.publish(th2_m)
        self.th3_pub.publish(th3_m)
        self.th4_pub.publish(th4_m)
        self.th5_pub.publish(th5_m)

    def inv_kinematics(self, Px, Py, Pz):
        a2, a3, d3, d4 = self.a[2], self.a[3], self.d[2], self.d[3]

        th1 = m.atan2(Py, Px) - m.atan2(d3, m.sqrt(Px**2 + Py**2 - d3**2))
        
        k = (Px**2 + Py**2 + Pz**2 - a2**2 - a3**2 - d3**2 - d4**2)/(2*a2)
        th3 = m.atan2(a3, d4) - m.atan2(k, m.sqrt(a3**2 + d4**2 - k**2))

        th2 = m.atan2((Px*m.cos(th1) + Py*m.sin(th1))*(a2*m.sin(th3) - d4) - Pz*(a3 + a2*m.cos(th3)),
            (a3 + a2*m.cos(th3))*(Px*m.cos(th1) + Py*m.sin(th1)) + (a2*m.sin(th3) - d4)*Pz) - th3
        
        r13 = self.Toriginefector.row(0).column(2)
        r23 = self.Toriginefector.row(1).column(2)
        r33 = self.Toriginefector.row(2).column(2)

        th4 = m.atan2()


if __name__ == "__main__":
    try:
        p560_c = puma560_control()
        rospy.spin()
    except (ROSInterruptException):
        sys.exit()