import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import math

def generate_atoms():
	r = 3.822  #r_min for argon
	x_position = []
	y_position = []
	atom_position = []
	coor = 0.0
	for i in range(20):
		for j in range(20):
			x_position.append(coor)
		y_position.append(coor)
		coor += r
	y_position *= 20
	atom_position.append(x_position)
	atom_position.append(y_position)

	return atom_position

def neighbor_list(atom_position):
	'''atom_position is the 
	(2,400) array that contains the positions of
	the 400 atoms block
	This function returns the list of neighbors map
	of each atom in the block
    And a numNeigh list that contains the number 
    of neighbor of each atom
	'''
	cut_off = 7.5
	neighbors = []
	numNeigh = []
	for i in range(400):
		lst = []
		for j in range(400):
			if j != i:
				r = math.sqrt((atom_position[0][i] - atom_position[0][j])**2 + (atom_position[1][i]-atom_position[1][j])**2)
				if r < cut_off:
					lst.append(j)
		neighbors.append(lst)
	for i in range(len(neighbors)):
		numNeigh.append(len(neighbors[i]))
	
	return neighbors, numNeigh
	
def plot_neighbors(atom_position):
	x = atom_position[0]
	y = atom_position[1]
	fig, ax = plt.subplots()
	atoms = plt.plot(x,y,'ro')
	plt.xlabel("x_position (Å)")
	plt.ylabel("y_position (Å)")
	circle = plt.Circle((atom_position[0][300],atom_position[1][30]),7.5,color = 'blue')
	ax.add_artist(circle)
	plt.show()

def calc_velocity(atom_position,neighbors):
	sigma = 3.405
	epsilon = 0.010323
	x = atom_position[0]
	xx = atom_position[0][:]
	y = atom_position[1]
	F_X = []
	F_Y = []

	kb = 8.617*(10**-5)
	dt = 10**-13
	vx = np.zeros(400)
	vy = np.zeros(400)
	velocity = np.zeros(400)
	KE = np.zeros(8000)
	m = 6.634*(10**-26)

	n = 0
	while n < 3500:
		#print("Iteration is " + str(n))
		F_X = []
		F_Y = []
		#autoco = 0
		for i in range(400):
			fx = 0
			fy = 0
			for j in range(len(neighbors[i])):
				r = math.sqrt((atom_position[0][neighbors[i][j]] - atom_position[0][i])**2 + (atom_position[1][neighbors[i][j]]-atom_position[1][i])**2)
				if x[neighbors[i][j]]-x[i] == 0:
					fx +=0
				else:
					fx += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((x[neighbors[i][j]]-x[i])/r)
				if y[neighbors[i][j]]-y[i] == 0:
					fy +=0
				else:
					fy += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((y[neighbors[i][j]]-y[i])/r)
			F_X.append(fx)
			F_Y.append(fy)
	
		for i in range(400):
			x[i] = x[i] + vx[i]*dt + dt**2/(2*m)*F_X[i]
			y[i] = y[i] + vy[i]*dt + dt**2/(2*m)*F_Y[i]

		neighbors, numNeigh = neighbor_list(atom_position)
		F_X1 = []
		F_Y1 = []

		for i in range(400):
			fx = 0
			fy = 0
			for j in range(len(neighbors[i])):
				r = math.sqrt((atom_position[0][neighbors[i][j]] - atom_position[0][i])**2 + (atom_position[1][neighbors[i][j]]-atom_position[1][i])**2)
				if x[neighbors[i][j]]-x[i] == 0:
					fx +=0
				else:
					fx += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((x[neighbors[i][j]]-x[i])/r)
				if y[neighbors[i][j]]-y[i] == 0:
					fy +=0
				else:
					fy += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((y[neighbors[i][j]]-y[i])/r)
			F_X1.append(fx)
			F_Y1.append(fy)

		K = 0
		for i in range(400):
			#vprev = math.sqrt((vx[i]**2+vy[i]**2))
			vx[i] = vx[i] + dt/(2*m)*(F_X1[i]+F_X[i])
			vy[i] = vy[i] + dt/(2*m)*(F_Y1[i]+F_Y[i])
			velocity[i] = math.sqrt((vx[i]**2+vy[i]**2))
			#autoco = autoco + vprev * velocity[i]
			K += 1/2*m*velocity[i]**2

		KE[n] = K

		n += 1

	n = 0
	while n < 4500:
		#print("Iteration is " + str(n))
		F_X = []
		F_Y = []
		for i in range(400):
			fx = 0
			fy = 0
			for j in range(len(neighbors[i])):
				r = math.sqrt((atom_position[0][neighbors[i][j]] - atom_position[0][i])**2 + (atom_position[1][neighbors[i][j]]-atom_position[1][i])**2)
				if x[neighbors[i][j]]-x[i] == 0:
					fx +=0
				else:
					fx += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((x[neighbors[i][j]]-x[i])/r)
				if y[neighbors[i][j]]-y[i] == 0:
					fy +=0
				else:
					fy += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((y[neighbors[i][j]]-y[i])/r)
			F_X.append(fx)
			F_Y.append(fy)
	
		for i in range(400):
			x[i] = x[i] + vx[i]*dt + dt**2/(2*m)*F_X[i]
			y[i] = y[i] + vy[i]*dt + dt**2/(2*m)*F_Y[i]

		neighbors, numNeigh = neighbor_list(atom_position)
		F_X1 = []
		F_Y1 = []
		PE = []

		for i in range(400):
			fx = 0
			fy = 0
			for j in range(len(neighbors[i])):
				r = math.sqrt((atom_position[0][neighbors[i][j]] - atom_position[0][i])**2 + (atom_position[1][neighbors[i][j]]-atom_position[1][i])**2)
				if x[neighbors[i][j]]-x[i] == 0:
					fx +=0
				else:
					fx += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((x[neighbors[i][j]]-x[i])/r)
				if y[neighbors[i][j]]-y[i] == 0:
					fy +=0
				else:
					fy += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((y[neighbors[i][j]]-y[i])/r)

			F_X1.append(fx)
			F_Y1.append(fy)

		K = 0
		for i in range(400):
			vx[i] = vx[i] + dt/(2*m)*(F_X1[i]+F_X[i])
			vy[i] = vy[i] + dt/(2*m)*(F_Y1[i]+F_Y[i])

			velocity[i] = math.sqrt((vx[i]**2+vy[i]**2))
			K += 1/2*m*velocity[i]**2
		KE[n+3500] = K
		if n % 500 == 0:
			vx = np.zeros(400)
			vy = np.zeros(400)
			velocity = np.zeros(400)

		n += 1

	#Voronoi polyhedra
	#volume = 3.822**2

	#hydro = []
	#P = 0
	#for i in range(400):
		#sum11 = 0
		#sum12 = 0
		#sum22 = 0
		#for j in range(len(neighbors[i])):
			#rx = atom_position[0][neighbors[i][j]] - atom_position[0][i]
			#ry = atom_position[1][neighbors[i][j]] - atom_position[1][i]
			#sum11 += F_X1[j]*rx
			#sum12 += F_X1[j]*ry
			#sum22 += F_Y1[j]*ry

		#sig11 = -1/volume*sum11
		#sig12 = -1/volume*sum12
		#sig22 = -1/volume*sum22
		#hy = -(sig11+sig22)/2
		#hydro.append(hy)
		#P = sum(hydro)
	#print("Total pressure of initial configuration equals to" + str(P))

	#x_l = np.arange(0,8000)
	#plt.plot(x_l,KE)
	#plt.xlabel("iterations")
	#plt.ylabel("Kinetic Energy")
	#plt.show()

	#fig = plt.scatter(x,y,c=hydro) # plot with coloring by the first column
	#a = plt.colorbar() # create colorbar
	#a.set_label("Hydrostatic Pressure (eV/${\AA^2}$)") # label colorbar
	#plt.title('Hydrostatic Pressure')
	#plt.xlabel('x (Å)')
	#plt.ylabel('y (Å)')
	#plt.show() # display the plot

	#calc_RDF(0.01,atom_position,neighbors)

	T = K/400/kb
	print("temperature" + str(T))

	Td = 5
	Total_P = []
	PE = np.zeros(4000)
	TotalE = np.zeros(4000)
	KE = np.zeros(4000)
	v0 = np.zeros(400)
	Total_P = []
	TE = []
	cvs = []

	while Td <= 80:
		cvv = []
		print("Td = " + str(Td))
		PE = np.zeros(4000)
		TotalE = np.zeros(4000)
		KE = np.zeros(4000)
		kappa = (Td/T)**0.5
		print("kappa is " + str(kappa))
		vx = vx * kappa
		vy = vy * kappa
		velocities = [[] for i in range(400)]
		n = 0
		pres = []
		#autoco = 0
		autocos = []
		while n < 4000:
			e = 0
			esquare = 0
			F_X = []
			F_Y = []
			for i in range(400):
				fx = 0
				fy = 0
				for j in range(len(neighbors[i])):
					r = math.sqrt((atom_position[0][neighbors[i][j]] - atom_position[0][i])**2 + (atom_position[1][neighbors[i][j]]-atom_position[1][i])**2)
					if x[neighbors[i][j]]-x[i] == 0:
						fx +=0
					else:
						fx += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((x[neighbors[i][j]]-x[i])/r)
					if y[neighbors[i][j]]-y[i] == 0:
						fy +=0
					else:
						fy += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((y[neighbors[i][j]]-y[i])/r)
				F_X.append(fx)
				F_Y.append(fy)
	
			for i in range(400):
				x[i] = x[i] + vx[i]*dt + dt**2/(2*m)*F_X[i]
				y[i] = y[i] + vy[i]*dt + dt**2/(2*m)*F_Y[i]
				
			neighbors, numNeigh = neighbor_list(atom_position)
			F_X1 = []
			F_Y1 = []
			for i in range(400):
				fx = 0
				fy = 0
				for j in range(len(neighbors[i])):
					r = math.sqrt((atom_position[0][neighbors[i][j]] - atom_position[0][i])**2 + (atom_position[1][neighbors[i][j]]-atom_position[1][i])**2)
					if x[neighbors[i][j]]-x[i] == 0:
						fx +=0
					else:
						fx += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((x[neighbors[i][j]]-x[i])/r)
					if y[neighbors[i][j]]-y[i] == 0:
						fy +=0
					else:
						fy += 24*epsilon/r*(-2*(sigma/r)**12+(sigma/r)**6)*((y[neighbors[i][j]]-y[i])/r)
					if r == 0:
						e += 0
					else:
						e += 4*epsilon*((sigma/r)**12-(sigma/r)**6)
						if n >= 2000:
							esquare += ((4*epsilon*((sigma/r)**12-(sigma/r)**6))/2)**2
				F_X1.append(fx)
				F_Y1.append(fy)
				#print(F_Y1)
			if n >= 2000:
				eavesquare = (e/400/2)**2
				esquare = esquare/400
				cv = (esquare-eavesquare)/(kb*Td**2)
				cvv.append(cv)
				


			PE[n]=e/2

			K = 0
			
			volume = 3.822**2

			hydro = []
			P = 0
			for i in range(400):

				#v0 = math.sqrt((vx[i]**2+vy[i]**2))

				vx[i] = vx[i] + dt/(2*m)*(F_X1[i]+F_X[i])
				vy[i] = vy[i] + dt/(2*m)*(F_Y1[i]+F_Y[i])
				velocity[i] = math.sqrt((vx[i]**2+vy[i]**2))
				if n >= 2000:
					velocities[i].append(velocity[i])

					sum11 = 0
					sum12 = 0
					sum22 = 0
					for j in range(len(neighbors[i])):
						rx = atom_position[0][neighbors[i][j]] - atom_position[0][i]
						ry = atom_position[1][neighbors[i][j]] - atom_position[1][i]
						sum11 += F_X1[j]*rx
						sum12 += F_X1[j]*ry
						sum22 += F_Y1[j]*ry

					sig11 = -1/volume*sum11
					sig12 = -1/volume*sum12
					sig22 = -1/volume*sum22
					hy = -(sig11+sig22)/2
					P += hy


				K += 1/2*m*velocity[i]**2

			if P != 0:
				pres.append(P)

			KE[n] = K
			TotalE[n] = KE[n]+PE[n]

			#print("iteration is " + str(n))
			if n % 500 == 0:
				T = K/400/kb
				kappa = (Td/T)**0.5
				vx = vx * kappa
				vy = vy * kappa
				
				#Voronoi polyhedra
				#volume = 3.822**2

				#hydro = []
				#P = 0
				#for i in range(400):
					#sum11 = 0
					#sum12 = 0
					#sum22 = 0
					#for j in range(len(neighbors[i])):
						#rx = atom_position[0][neighbors[i][j]] - atom_position[0][i]
						#ry = atom_position[1][neighbors[i][j]] - atom_position[1][i]
						#sum11 += F_X1[j]*rx
						#sum12 += F_X1[j]*ry
						#sum22 += F_Y1[j]*ry

					#sig11 = -1/volume*sum11
					#sig12 = -1/volume*sum12
					#sig22 = -1/volume*sum22
					#hy = -(sig11+sig22)/2
					#hydro.append(hy)
					#P = sum(hydro)
					#pres.append(P)

			n += 1

		
		#print("autoco calculating")
		#for i in range(1000):
			#autoco = 0
			#for j in range(1000):
				#for k in range(400):
					#v0 = velocities[k][j]
					#vt = velocities[k][j+i]
					#autoco = autoco + v0*vt

			#autoco_new = autoco/(400*1000)
			#autocos.append(autoco_new)
			#print(i)

		#Voronoi polyhedra
		#volume = 3.822**2

		#hydro = []
		#P = 0
		#for i in range(400):
			#sum11 = 0
			#sum12 = 0
			#sum22 = 0
			#for j in range(len(neighbors[i])):
				#rx = atom_position[0][neighbors[i][j]] - atom_position[0][i]
				#ry = atom_position[1][neighbors[i][j]] - atom_position[1][i]
				#sum11 += F_X1[j]*rx
				#sum12 += F_X1[j]*ry
				#sum22 += F_Y1[j]*ry

			#sig11 = -1/volume*sum11
			#sig12 = -1/volume*sum12
			#sig22 = -1/volume*sum22
			#hy = -(sig11+sig22)/2
			#hydro.append(hy)
			#P = sum(hydro)

		cv = sum(cvv)/len(cvv)
		cvs.append(cv)

		P = sum(pres)/len(pres)
		print("Total pressure equals to" + str(P))
		Total_P.append(P)

		TE.append(TotalE[3999])
		print('total energy is '+ str(TotalE[3999]))

		#x_l = np.arange(0,1000)
		#plt.plot(x_l, autocos)
		#plt.xlabel("Time (10E-13 s)")
		#plt.ylabel("Autocorrelation function (${\AA^2}$/fs^2)")
		#plt.show()

		#x_l = np.arange(0,4000)
		#plt.plot(x_l,KE)
		#plt.xlabel("iterations")
		#plt.ylabel("Kinetic Energy (eV)")
		#plt.show()

		#fig = plt.scatter(x,y,c=hydro) # plot with coloring by the first column
		#a = plt.colorbar() # create colorbar
		#a.set_label("Hydrostatic Pressure (eV/${\AA^2}$)") # label colorbar
		#plt.title('Hydrostatic Pressure')
		#plt.xlabel('x (Å)')
		#plt.ylabel('y (Å)')
		#plt.show() # display the plot

		#fig = plt.scatter(x,y,c=hydro) # plot with coloring by the first column
		#a = plt.colorbar() # create colorbar
		#a.set_label("Hydrostatic Pressure (eV/${\AA^2}$)") # label colorbar
		#plt.title('Hydrostatic Pressure')
		#axes = plt.gca()
		#axes.set_xlim([0,80])
		#axes.set_ylim([0,80])
		#plt.xlabel('x (Å)')
		#plt.ylabel('y (Å)')
		#plt.show() # display the plot

		#calc_RDF(0.01,atom_position,neighbors)

		Td += 5

	temperatures = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]
	print("Pressures are " )
	print(Total_P)
	plt.plot(temperatures,Total_P)
	plt.xlabel("Temperature (K)")
	plt.ylabel("Total Pressure (eV/${\AA^2}$)")
	plt.show()

	#print("Energies are ")
	#print(TE)
	#plt.plot(temperatures,TE)
	#plt.xlabel("Temperature (K)")
	#plt.ylabel("Total Energy (eV)")
	#plt.show()

	print("Cvs are ")
	print(cvs)
	plt.plot(temperatures,cvs)
	plt.xlabel("Temperature (K)")
	plt.ylabel("Heat Capacity (eV/K)")
	plt.show()

	return atom_position, neighbors
	
def calc_RDF(delta_r,atom_position,neighbors):
	V = 72.618**2
	N = 400
	pho = N/V
	pi = 3.1415926
	r = 0.5
	gr = []
	lst = []
	while r < 15:
		count = 0
		nr = 0
		for i in range(400):
			for j in range(400):
				radius = math.sqrt((atom_position[0][i] - atom_position[0][j])**2 + (atom_position[1][i]-atom_position[1][j])**2)
				if r < radius < r + delta_r:
					count += 1
					nr += 1/N*count
		gr.append(1/pho*nr/(2*pi*r*delta_r))
		lst.append(r)
		r += 0.01

	plt.plot(lst,gr)
	plt.title("RDF vs Radius for delta_r = %s" %(delta_r))
	plt.xlabel("Radius (Å)")
	plt.ylabel("RDF g(r)")
	plt.show()




def main():
	atom_position = generate_atoms()
	neighbors, numNeigh = neighbor_list(atom_position)
	atom_position1, neighbors1 = calc_velocity(atom_position,neighbors)
	#calc_RDF(0.01,atom_position1,neighbors1)

if __name__ == '__main__':
	main()