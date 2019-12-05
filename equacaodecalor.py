# ***************************************************
# Transfer of Heat and Mass homework 1
# author Caio E C Oliveira
# date October 30th, 2019.
# Python version 3.6.8
# Description:
	# This script uses the finite differeces method to solve the Heat Conduction Equation.
	# The temperature is represented by a list of lists. It will be refered to as matrix.
	# The matrix has its rows corresponding to a snapshot of the bar in a certain time, 
	# whereas its columns represent each part of it, separated into N pieces.
	# 	  	x0 		x1 		x2 ...
	# t0:	T(0,0)	T(0,1)	T(0,2)
	# t1: 	T(1,0)	T(1,1)	T(1,2)
	# ...

	# Von Neuman Stability Analysis
	# Also, note that alpha must be
	# 	:= deltat/(deltax^2) <= 0.5 
# ***************************************************
# Modules import
import numpy as np
import matplotlib.pyplot as plt

# Functions ***************************************************
def numsolver(T, alpha, ktotal, N):
	for k in range (0, ktotal - 1):
		for i in range(1, N - 1): # N - 1 because we already know the value of temperature at i = N
			T[k+1,i] = T[k,i] + alpha*(T[k, i + 1] - 2*T[k,i] + T[k, i-1]) 
	return T

def showresults(T, x_axis, time_stamps):
	for i in range(len(time_stamps)):
		plt.plot(x_axis, T[time_stamps[i]])
	plt.xlabel('Bar length (m)') 
	plt.ylabel('Bar temperature (Celsius)') 
	plt.title('Temperature vs Bar length') 
	plt.legend()
	plt.show()

def T_analitica(x,t):
	return np.exp(-np.pi*np.pi*t)*np.sin(np.pi*x/2)

def problem_slot(L, deltat, ttotal, N, PROBLEM_OPTION):
	# Defining constants
	ktotal = 1 + int(ttotal/(deltat))
	deltax = L/N
	x_axis = np.arange(0, L, deltax)
	alpha  = deltat/(deltax**2)
	T 	   = np.zeros ([ktotal, N])

	if (PROBLEM_OPTION == 1):
		print ('Problem 1 initial conditions.')
		try:
			input('Press any key to continue.\n')
		except:
			pass
		T[0, 1:] 	= 1 	# 0 < x < 1, t = 0
		T[:, 0] 	= 0 	# x = 0, t >= 0
		T[:, -1]	= 0 	# x = 1, t >= 0
	elif (PROBLEM_OPTION == 2):
		print ('Problem 2 initial conditions.')
		try:
			input('Press any key to continue.\n')
		except:
			pass
		T[0, 1:] 	= 0 	# 0 < x < 1, t = 0
		T[:, 0] 	= 1 	# x = 0, t >= 0
		T[:, -1]	= 0 	# x = 1, t >= 0
	elif (PROBLEM_OPTION == 3):
		print ('Problem 3 initial conditions.')
		try:
			input('Press any key to continue.')
		except:
			pass
		T[0, 1:] 	= np.sin((np.pi/2)*x_axis[1:]) 	# 0 < x < 2, t = 0
		T[:, 0] 	= 0 	# x = 0, t >= 0
		T[:, -1] 	= 0 	# x = 1, t >= 0
		for t in range (len(x_axis)):
			plt.plot(x_axis, T_analitica(x_axis,t))
		plt.xlabel('Bar length (m)') 
		plt.ylabel('Bar temperature (Celsius)') 
		plt.title('Temperature vs Bar length analytical solution') 
		plt.legend()
		plt.show()
	else:
		print('\n\nOption not found...')
		exit()
	print ('Deltat =', deltat)
	print ('Deltax =', deltax)
	print ('Alpha =', alpha)
	print ('Bar length(L) =', L, '\nNumber of divisions(N) =', N)
	print ('Deltat =', deltat, '\nDeltax = L/N =', deltax)
	print ('Total simulation run time =', ttotal, '\n*********')

	T = numsolver(T = T, alpha = alpha, ktotal = ktotal, N = N)
	time_stamps = [0, int(ttotal/3), int(2*ttotal/3), -1]
	showresults(T = T, x_axis = x_axis, time_stamps = time_stamps)

# Main
problem_slot(L = 1, deltat = 5e-3, ttotal = 100, N = 20, PROBLEM_OPTION = 1)
problem_slot(L = 1, deltat = 5e-3, ttotal = 100, N = 20, PROBLEM_OPTION = 2)
problem_slot(L = 2, deltat = 1e-4,  ttotal = 100, N = 100, PROBLEM_OPTION = 3)
