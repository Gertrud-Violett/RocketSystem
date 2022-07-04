# -*- coding: utf-8 -*-
# ======
# 複合材の剛性、強度の計算
# 任意のレイアップの各配置角における剛性強度の計算を行なう。
#
# Simple code for calulating composite material stiffness and strength
# Classic Laminate Theory and Tsai-Wu Failure Criterion
#
# 使い方：def layupに各物性値を直打ちしてお使いください
# Usage: Enter variables in def layup
#
# Copyright (c) 2019 K.Trude
# This code is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
# ======

#標準セット
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')


#read input file

case_name = 'test' #enter case name
if not os.path.exists(case_name):
    os.mkdir(case_name)

layup_origin = np.array([0,0,0,0,0,0,0,0])    #積層構成の読込 Orientation of the fibres for each layer/ply of the laminate
num_layup = len(layup_origin)

def layup(orientation): #物性値の読込
    return layup_origin + np.full(num_layup, orientation)

#積層理論剛性パラメータ
E1 = np.full(num_layup, 125.)       # Axial stiffness of each ply [GPa]
#E1 = np.array([125,125,125])                 # (Use this instead if layer properties are not the same for all layers)

E2 = np.full(num_layup, 9.)         # transverse stiffness of each ply [GPa]
#E2 = np.array([16.2,16.2.16,2])                 # (Use this instead if layer properties are not the same for all layers)

v12 = np.full(num_layup, 0.33)      # Poisson's ratio
#v12 = np.array([0.278,0.278,0.278])             # (Use this instead if layer properties are not the same for all layers)
v21 = (v12*E2)/E1                                # (Since the compliance matrix is symmetric, v12/E1=v21/E2)

E45 = np.full(num_layup, 12.)       # (Use for Manual entry to calculate G12)
nu_45 = np.full(num_layup, 0.31)    # (Use for Manual entry to calculate G12)
G12_man = E45/2/(1+nu_45)

G12 = np.full(num_layup, 4.5)       # Shear modulus [GPa]
#G12 = np.array([5.83,5.83,5.83])                # (Use this instead if layer properties are not the same for all layers)
print('G12 calculated value =' + str(G12_man[0]))

h0 = 0.2
h = np.full(num_layup, h0)
# h = np.array([h0,h0,h0,h0,h0,h0,h0,h0])                   #height of each ply [mm] 枚数を合わせる


#強度則パラメータ
F_Lt = np.full(num_layup, 2700)       # 0 deg Tensile Strength of each ply
#F_Lt = np.array([2000,2000,2000]) 

F_Lc = np.full(num_layup, 680)       # 0 deg Compressive Strength of each ply
#F_Lc = np.array([2000,2000,2000]) 

F_Tt = np.full(num_layup, 57.)       # 90 deg Tensile Strength of each ply
#F_Tt = np.array([2000,2000,2000]) 

F_Tc = np.full(num_layup, 130.)       # 90 deg Compressive Strength of each ply
#F_Tc = np.array([2000,2000,2000]) 

F_LTs = np.full(num_layup, 88.)       # Shear Strength of each ply (45°引張試験の実測値)
#F_LTs = np.array([2000,2000,2000]) 

# Looping through each layer of the laminate
def calc_Q(orientation):
    # Initiating Q as a list that can contain the stiffness matrices transformed into the global coordinates for each layer
    Q=[]

    for i in range(num_layup):
        theta = layup(orientation)[i]*np.pi/180 # Current ply angle changed into radians
        # assigning the current material properties to temporary variables (could also be used directly):
        E1l=E1[i]
        E2l=E2[i]
        v12l=v12[i]
        v21l=v21[i]
        G12l=G12[i]
        #print('Ply no. ' + str(i+1) + ': theta='+ str(layup(orientation)[i]) +', E1l=' + str(E1l), ' E2l=' + str(E2l))
        # Establishing current local stiffness matrix, Ql:
        Ql = 1/(1-v12l*v21l)*np.array([[E1l,v21l*E1l,0],\
                                       [v12l*E2l,E2l,0],\
                                       [0,0,G12l*(1-v12l*v21l)]])
        # Transformation matrix:
        T=np.array([[np.cos(theta)**2,np.sin(theta)**2,-2*np.sin(theta)*np.cos(theta)], \
                    [np.sin(theta)**2,np.cos(theta)**2,2*np.sin(theta)*np.cos(theta)],\
                    [np.sin(theta)*np.cos(theta),-np.sin(theta)*np.cos(theta),np.cos(theta)**2-np.sin(theta)**2]])
        # Adding the current stiffness matrix in the global coordinate system to the Q-list variable:
        Q.append(np.dot(np.dot(T,Ql),np.transpose(T)))
        #print(Q[i]) #debug function
    return Q
    
Q = calc_Q(0.)

def calc_A(Q, h, num_layup):
    A=np.zeros((3,3))   #specify layup ply no.
    for i in range(num_layup):
        A=A+Q[i]*h[i]
    return A
A = calc_A(Q, h, num_layup)

def calc_Qstar_and_Sstar(A, h): #剛性行列の計算def
    Qstar = A/sum(h)
    Sstar = np.linalg.inv(Qstar)
    Ex = 1/Sstar[0,0]
    Ey = 1/Sstar[1,1]
    Gxy = 1/Sstar[2,2]
    vxy=-Sstar[1,0]/Sstar[0,0]
    Gxy_check = Ex/2/(1+vxy)
    return Qstar, Sstar, Ex, Ey, Gxy, vxy, Gxy_check


Qstar, Sstar, Ex, Ey, Gxy, vxy, Gxy_check = calc_Qstar_and_Sstar(A, h)

# Printing out the matrices
print('Q* Matrix')
print(Qstar)
print('S* Matrix')
print(Sstar)
print('Stiffness Ex='+ str(Ex))
print('Stiffness Ey='+ str(Ey))
print('shear modulus Gxy='+ str(Gxy))
print('Stiffness Gxy(等方材のみ、検算用)='+ str(Gxy_check))
print('poisson ratio vxy='+str(vxy))
Ex0 = Ex
Ey0 = Ey
vxy0 = vxy

#Parameters for stress matrix F
F12_star = -0.5 #相互干渉項。よくわからなかったら-0.5にしておく
F11 = 1/(F_Lt*F_Lc)
F22 = 1/(F_Tt*F_Tc)
F12 = F12_star*np.sqrt(F11*F22)
F66 = 1/F_LTs**2
F1 = 1/F_Lt-1/F_Lc
F2 = 1/F_Tt-1/F_Tc

print('F11='+str(F11))
print('F22='+str(F22))
print('F12='+str(F12))
print('F66='+str(F66))
print('F1='+str(F1))
print('F2='+str(F2))
print('')

def calc_sigma(orientation): #破断応力の計算def (Tsai-Wu Failure Criterion破壊則)
    # Initiating list that can contain the stress matrices transformed into the global coordinates for each layer
    FA=[]
    FB=[]
    sigma_tl=[]
    sigma_cl=[]
    # Looping through each layer of the laminate
    for i in range(num_layup):
        theta = layup(orientation)[i]*np.pi/180 # Current ply angle changed into radians
        # assigning the current material properties to temporary variables (could also be used directly):
        F11l=F11[i]
        F22l=F22[i]
        F12l=F12[i]
        F66l=F66[i]
        F1l=F1[i]
        F2l=F2[i]
        #print('Ply no. ' + str(i+1) + ': theta='+ str(layup(orientation)[i]))
       
        FB = F11l*np.cos(theta)**4+F22l*np.sin(theta)**4+(2*F12l+F66l)*(np.cos(theta)**2*np.sin(theta)**2)
        FA = F1l*np.cos(theta)**2+F2l*np.sin(theta)**2
        sigma_x_t = (-FA+np.sqrt(FA**2+4*FB))/2/FB
        sigma_x_c = (FA+np.sqrt(FA**2+4*FB))/2/FB
        sigma_tl.append(sigma_x_t)
        sigma_cl.append(sigma_x_c)
        #print('Tensile Strength = %d MPa' %(sigma_x_t))
        #print('Compressive Strength = %d MPa'  %(sigma_x_c))
        #print('')

    return sigma_x_t, sigma_x_c, sigma_tl, sigma_cl

sigma_x_t, sigma_x_c, sigma_tl, sigma_cl = calc_sigma(0.)

h0 = 0.2
h = np.full(num_layup, h0)
# h = np.array([h0,h0,h0,h0,h0,h0,h0,h0])                   #height of each ply [mm] 枚数を合わせる
thickness = np.sum(h)
f = open(case_name + '/output_' + case_name + '.txt', 'w')
sys.stdout = f

print('Layup pattern'  + str(layup_origin))
print('Layup Total Thickness = %.2f mm'  %(thickness))
print('Combined Tensile Strength Ft (0 deg) = %d MPa' %(sum(sigma_tl*h)/thickness))
print('Combined Compressive Strength Fc (0 deg)= %d MPa'  %(sum(sigma_cl*h)/thickness))
print('Combined 0 deg Elastic Modulus Ex (0 deg) = %.1f'  %(Ex0))
print('Combined 90 deg Elastic Modulus Ey (0 deg) = %.1f'  %(Ey0))
print('poisson ratio (0 deg) v = %.3f'  %(vxy0))


# グラフ出力したい値を用意する
Ex_list = []
Ey_list = []
Total_Tensile_Strength_list = []
Total_Compressive_Strength_list = []

for orientation in range(91):
    Q = calc_Q(orientation)
    A = calc_A(Q, h, num_layup)
    Qstar, Sstar, Ex, Ey, Gxy, vxy, Gxy_check = calc_Qstar_and_Sstar(A, h)
    sigma_x_t, sigma_x_c, sigma_tl, sigma_cl = calc_sigma(orientation)
    Ex_list.append(Ex)
    Ey_list.append(Ey)
    Total_Tensile_Strength_list.append(sum(sigma_tl*h)/thickness)
    Total_Compressive_Strength_list.append(sum(sigma_cl*h)/thickness)

    orientation = range(91)

plt.figure()
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.xlabel('Orientation [deg]')
plt.ylabel('Stiffness [GPa]')
plt.plot(orientation, Ex_list)
plt.title("Orientaion - Layup Stiffness Ex")
plt.savefig(case_name + "/Stiffness_Ex.png")

plt.figure()
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.xlabel('Orientation [deg]')
plt.ylabel('Stiffness [GPa]')
plt.plot(orientation, Ey_list)
plt.title("Orientaion - Layup Stiffness Ey")
plt.savefig(case_name + "/Stiffness_Ey.png")

plt.figure()
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.xlabel('Orientation [deg]')
plt.ylabel('Strength [MPa]]')
plt.plot(orientation, Total_Tensile_Strength_list)
plt.title("Orientaion - Layup Tensile Strength")
plt.savefig(case_name + "/Strength_Tensile.png")

plt.figure()
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.xlabel('Orientation [deg]')
plt.ylabel('Strength [MPa]]')
plt.plot(orientation, Total_Compressive_Strength_list)
plt.title("Orientaion - Layup Compressive Strength")
plt.savefig(case_name + "/Strength_Compressive.png")