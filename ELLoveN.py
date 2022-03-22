###########################################
##  -----------------------------------  ##
##  | ELLoveN,  an Elastic Load Love  |  ##
##  | Numbers Calculator for a SNREI  |  ##
##  | incompressible Earth model.     |  ##
##  | G. Del Piccolo                  |  ##
##  -----------------------------------  ##
###########################################

# --- Import libraries ---------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import sys

# -----------------------------------------------------------------------

# --- Defining functions for the ODEs -----------------------------------

# --- Core propagation matrix
def Cprop(r,rho,g,n):
	m1 = np.array([[(4*np.pi*G*rho/g) -(n+1)/r, 1],
               [8*np.pi*G*rho*(n-1)/(g*r), ((n-1)/r)-(4*np.pi*G*rho/g)]])
	return m1

# --- Mantle propagation matrix
def Mprop(r,rho,mu,g,a,n):
	m2 = np.array([[-2/r,n*(n+1)/r,0,   0,0,0],
               [-1/r,      1/r,0,1/mu,0,0],
 [((3*mu/r)-rho*g)*4/r,((6*mu/r)-rho*g)*(-n*(n+1)/r),0,n*(n+1)/r,-rho*(n+1)/r,rho],
 [(rho*g-6*mu/r)/r,-(1/(r*r))*(2*mu-n*(n+1)*4*mu),-1/r,-3/r,rho/r,0],
 [-4*np.pi*(G)*rho,0,0,0,-(n+1)/r,1],
 [-4*np.pi*(G)*rho*(n+1)/r,4*np.pi*(G)*rho*n*(n+1)/r,0,0,0,(n-1)/r]])
	return m2

# Runge-Kutta2 integration algorithm for core and mantle

def RK2_core(y,steps,size,radius,dr,n,rho):

	for i in range(steps):
		prop = np.zeros((size,size),dtype=np.float64)
		gr = (4/3)*np.pi*G*rho*radius[i]
		prop = Cprop(radius_core[i],rho_c,gr,n)
		y_ct = np.zeros(size,dtype=np.float64)

		for j in range(size):
			prod = 0 
			for k in range(size):
				prod += prop.item(j,k)*y[(i)*size+(k)]
			y_ct[(j)] = dr*prod/2 + y[(i)*size+(j)] 
			gr = (4/3)*np.pi*G*rho_c*(radius[i]+dr_c/2)
		prop = Cprop(radius[i]+dr/2,rho,gr,n)

		for j in range(size):
			prod = 0 
			for k in range(size):
				prod += prop.item(j,k)*y_ct[(k)]
			y[(i+1)*size+(j)] = dr*prod + y[(i)*size+(j)] 

	return y


def RK2_mantle(y1,y2,y3,steps,size,radius,dr,n,rho,mu,g):

	for i in range(steps):
		prop = np.zeros((size,size),dtype=np.float64)
		prop = Mprop(radius[i],rho[i],mu[i],g[i],a,n)
		y1t = np.zeros(size,dtype=np.float64)
		y2t = np.zeros(size,dtype=np.float64)
		y3t = np.zeros(size,dtype=np.float64)
		for j in range(size):
			prod1 = 0 
			prod2 = 0
			prod3 = 0
			for k in range(size):
				prod1 += prop.item(j,k)*y1[(i)*size+(k)]
				prod2 += prop.item(j,k)*y2[(i)*size+(k)]
				prod3 += prop.item(j,k)*y3[(i)*size+(k)]
			y1t[(j)] = dr*prod1/2 + y1[(i)*size+(j)]
			y2t[(j)] = dr*prod2/2 + y2[(i)*size+(j)]
			y3t[(j)] = dr*prod3/2 + y3[(i)*size+(j)] 

			prop = Mprop(radius[i]+dr/2,rho[i],mu[i],g[i],a,n)
			for j in range(size):
				prod1 = 0 
				prod2 = 0
				prod3 = 0
				for k in range(size):
					prod1 += prop.item(j,k)*y1t[(k)]
					prod2 += prop.item(j,k)*y2t[(k)]
					prod3 += prop.item(j,k)*y3t[(k)]
				y1[(i+1)*size+(j)] = dr*prod1 + y1[(i)*size+(j)]
				y2[(i+1)*size+(j)] = dr*prod2 + y2[(i)*size+(j)]
				y3[(i+1)*size+(j)] = dr*prod3 + y3[(i)*size+(j)]

	y_lin = np.zeros((3,size*(steps+1)),dtype=np.float64)
	y_lin[0,:] = y1
	y_lin[1,:] = y2
	y_lin[2,:] = y3

	return y_lin

# -----------------------------------------------------------------------

# --- Declaring variables -----------------------------------------------

max_deg = int(sys.argv[2]) + 1    # Max harmonic degree 
min_deg = int(sys.argv[1])        # Min harmonic degree

# --- Scales for mass, time and length

Ma = 4/(3)*np.pi*np.power(6.371*np.power(10,6),3)*5517
T = np.sqrt(np.power(6.371*np.power(10,6),3)/(Ma*(6.67*(10**(-11)))))
L = 6.371*np.power(10,6)

# --- Non-dimensionalized quantities 

Me = 1                             # Earth Mass 
rho_c = 10750*np.power(L,3)/Ma     # Core's density
g0 = 9.82 * np.power(T,2) /(L)     # Gravity field at surface
a = 1                              # Earth's radius
b = 0.3480/0.6371                  # Core's radius
mrad = 0.34/0.6371                 # Initialization radius inside core
G = 1                              # Gravitational constant
steps_core = 50                    # Iterations inside Core
steps_mantle = 200                 # Iterations inside Mantle + Lith.

# --- Elastic Load Love Numbers arrays 

h_n = np.zeros(max_deg-min_deg)    # Vertical displacement
l_n = np.zeros(max_deg-min_deg)    # Horizontal displacement
k_n = np.zeros(max_deg-min_deg)    # Geoid displacement
deg = np.zeros(max_deg-min_deg)  # Harmonic degrees

# --- Output files

f1 = open("h_n.dat","w")
f2 = open("l_n.dat","w")
f3 = open("k_n.dat","w") 

# -----------------------------------------------------------------------


# --- Initialize radius arrays inside Core and Mantle+Lith --------------

radius_core = np.linspace(mrad,b,steps_core)
dr_c = (b-mrad)/(steps_core)  # Core Spacing 

radius_mantle = np.linspace(b,a,steps_mantle)
dr_m = (a-b)/(steps_mantle)   # Mantle Spacing 


#------------------------------------------------------------------------
print(" ")
print("|--------------------------------------|")
print("|--------------------------------------|")
print("||   ELLoveN,  an Elastic Load Love   ||")
print("||   Numbers Calculator for a SNREI   ||")
print("||   incompressible Earth.            ||")
print("||                                    ||")
print("||                    G. Del Piccolo  ||")
print("|--------------------------------------|")
print("|--------------------------------------|")
print(" ")

print("Loading Earth Model...")
ln = 0
f = open('Earth_Model.txt', 'r')
minrad, maxrad, rigidity, density, gravity = [], [], [], [], []
for line in f.readlines():
	ln += 1
	fields = line.split(',')
	minrad.append(float(fields[0]))
	maxrad.append(float(fields[1]))
	rigidity.append(float(fields[2]))
	density.append(float(fields[3]))
	gravity.append(float(fields[4]))
f.close()

# This parameter determines the profile of rigidity, density and gravity inside the mantle. Set it to False in order to have shells with uniform properties inside, set it to True in order to have a linear profile for the properties in the mantle. 
linear = sys.argv[3]
mu = np.zeros(steps_mantle)
rho_m = np.zeros(steps_mantle)
g = np.zeros(steps_mantle)

if(linear == "constant"):
	for i in range(steps_mantle):
		x = radius_mantle[i]
		for j in range(1,ln,1):
			if ((x>=(minrad[j]/6371))&(x<=(maxrad[j]/6371))): 
				mu[i] = rigidity[j]*(10**11) * L * np.power(T,2) / (Ma)  
				rho_m[i] = density[j]*np.power(L,3)/Ma
				g[i] = gravity[j] * np.power(T,2) /(L) 

if(linear == "linear"): 
	for i in range(steps_mantle):
		x = radius_mantle[i]
		for j in range(1,ln,1):
			if ((x>=(minrad[j]/6371))&(x<=(maxrad[j]/6371))&(j<(ln-1))): 
				mu[i] = (rigidity[j]+((rigidity[j+1]-rigidity[j])*(x-(minrad[j]/6371)))/(maxrad[j]/6371-minrad[j]/6371))*(10**11)*L*np.power(T,2)/Ma 
				rho_m[i] = (density[j]+((density[j+1]-density[j])*(x-(minrad[j]/6371)))/(maxrad[j]/6371-minrad[j]/6371))*np.power(L,3)/Ma
				g[i] = (gravity[j]+((gravity[j+1]-gravity[j])*(x-(minrad[j]/6371)))/(maxrad[j]/6371-minrad[j]/6371))*np.power(T,2)/L
			if ((x>=(minrad[j]/6371))&(x<=(maxrad[j]/6371))&(j==(ln-1))):
				mu[i] = rigidity[j]*(10**11) * L * np.power(T,2) / (Ma)  
				rho_m[i] = density[j]*np.power(L,3)/Ma
				g[i] = gravity[j] * np.power(T,2) /L

print("Earth LOADED")
print(" ")

# --- Harmonic Loop -----------------------------------------------------

for n in range (min_deg,max_deg,1):

	print("Processing: " + str(int((n-min_deg)/(max_deg-min_deg-1)*100))+"%",end="\r")

# --- Initialization harmonic components of the fields 
	U_n = 0
	V_n = 0
	P_n = 0

# --- Initialization independent solutions outside core 
	y_lin = np.zeros((3,6*(steps_mantle+1)),dtype=np.float64)
	y1 = np.zeros(6*(steps_mantle+1),dtype=np.float64)
	y2 = np.zeros(6*(steps_mantle+1),dtype=np.float64)
	y3 = np.zeros(6*(steps_mantle+1),dtype=np.float64)

# --- Initialization of the solution inside the core
	y_c = np.zeros(2*(steps_core+1),dtype=np.float64)
	y_c[0]  = mrad**n
	y_c[1]  = 2*(n-1)*mrad**(n-1)

# --- Propagation solution inside the core using Runge-Kutta 2 method
	y_c = RK2_core(y_c,steps_core,2,radius_core,dr_c,n,rho_c)

# --- Initialization the solution at the CMB
	y1[0] = -y_c[steps_core*2]/(4*np.pi*G*rho_c*b/3)
	y1[1] = 0
	y1[2] = 0
	y1[3] = 0
	y1[4] = y_c[steps_core*2+0]
	y1[5] = y_c[steps_core*2+1]

	y2[0] = 0
	y2[1] = 1
	y2[2] = 0
	y2[3] = 0
	y2[4] = 0 
	y2[5] = 0

	y3[0] = 1
	y3[1] = 0
	y3[2] = rho_c*4*np.pi*G*rho_c*b/3
	y3[3] = 0
	y3[4] = 0
	y3[5] = 4*np.pi*G*rho_c


# --- Propagation solution outside the core using Runge-Kutta 2 method
	y_lin = RK2_mantle(y1,y2,y3,steps_mantle,6,radius_mantle,dr_m,n,rho_m,mu,g)
	y1 = y_lin[0,:]
	y2 = y_lin[1,:]
	y3 = y_lin[2,:]
	

# --- Boundary conditions at the surface --------------------------------

	M = np.array([[y1[6*(steps_mantle+1)-1-3],y2[6*(steps_mantle+1)-1-3],y3[6*(steps_mantle+1)-1-3]],
              [y1[6*(steps_mantle+1)-1-2],y2[6*(steps_mantle+1)-1-2],y3[6*(steps_mantle+1)-1-2]],
              [y1[6*(steps_mantle+1)-1],y2[6*(steps_mantle+1)-1],y3[6*(steps_mantle+1)-1]]])
	invM = np.zeros((3,3),dtype=np.float64)
	invM = np.linalg.inv(M)
	B = np.array([-g0*(2*n+1)/(4*np.pi*a*a),0,-g0*(2*n+1)/Me])
	c = np.zeros(3,dtype=np.float64)
	for p in range(3):
		c[p] = 0
		for q in range(3):
			c[p] += invM.item(p,q)*B[q]
	U_n = c[0]*y1[6*(steps_mantle+1)-1-5] + c[1]*y2[6*(steps_mantle+1)-1-5] + c[2]*y3[6*(steps_mantle+1)-1-5]
	V_n = c[0]*y1[6*(steps_mantle+1)-1-4] + c[1]*y2[6*(steps_mantle+1)-1-4] + c[2]*y3[6*(steps_mantle+1)-1-4]
	P_n = c[0]*y1[6*(steps_mantle+1)-1-1] + c[1]*y2[6*(steps_mantle+1)-1-1] + c[2]*y3[6*(steps_mantle+1)-1-1]
	h_n[n-min_deg] = (U_n * Me/a)
	l_n[n-min_deg] = (V_n * Me/a)
	k_n[n-min_deg] = (-1 - P_n * Me / a / g0)
	deg[n-min_deg] = n

# --- Write to output file

	f1 = open("h_n.dat","a")
	f2 = open("l_n.dat","a")
	f3 = open("k_n.dat","a") 
	ln_h = "{},{} \n".format(int(deg[n-min_deg]),h_n[n-min_deg])
	ln_l = "{},{} \n".format(int(deg[n-min_deg]),l_n[n-min_deg])
	ln_k = "{},{} \n".format(int(deg[n-min_deg]),k_n[n-min_deg])
	f1.writelines(ln_h)
	f2.writelines(ln_l)
	f3.writelines(ln_k)

f1.close()
f2.close()
f3.close()

# -----------------------------------------------------------------------

# --- Plot --------------------------------------------------------------

fig0,ax0 = plt.subplots(2,2,figsize=(8,6))
ax0[0,0].plot(radius_mantle*L/(10**6),mu/(L * np.power(T,2)*(10**9) / (Ma)),color="red")
ax0[0,0].set_title("rigidity (GPa) vs radius (km) in Mantle+Lith.",fontsize=8)
ax0[0,1].plot(radius_mantle*L/(10**6),rho_m/(np.power(L,3)/Ma),color="red")
ax0[0,1].set_title("density ($kg$ $m^{-3}$) vs radius (km) in Mantle+Lith.",fontsize=8)
ax0[1,0].plot(radius_mantle*L/(10**6),g/(np.power(T,2)/L),color="red")
ax0[1,0].set_title("gravity ($m$ $s^{-2}$) vs radius (km) in Mantle+Lith.",fontsize=8)

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(deg,h_n,"ok",markersize=3,label="h_n")
ax.plot(deg,l_n,"om",markersize=3,label="l_n")
ax.plot(deg,k_n,"ob",markersize=3,label="k_n")
ax.set_xscale("log")
ax.set_title("Incompressible Elastic Earth with Inviscid Core")
ax.set_xlabel("n",fontsize=15)
ax.tick_params(axis="x",labelsize=14)
ax.tick_params(axis="y",labelsize=14)
ax.legend(loc = "lower left")
ax.set_xlabel("n")
ax.grid(linestyle="--",color = "gray")

# -----------------------------------------------------------------------

plt.show()



