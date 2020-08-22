"""
Front tracking code for conservation laws with discontinuous flux

There are two numerical experiments being performed:
Experiment 1: Initial datum is a Riemann initial datum
			  Flux switches from Transport to Burgers across x=0

Experiment 2: Initial datum is a smooth bump
			  Flux switches from Burgers to Transport across x=0

To get a table with the convergence rates for each experiment uncomment the corresponding line in main
"""

import numpy as np

import scipy.integrate as integrate

from scipy.optimize import brentq

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import lines


class Method(object):
	""" Class that sets up the parameters needed for Front Tracking and the solution method itself """
	def __init__(self,u0,a,b,T,delta,M,g,f,finverse):
		self.u0,self.a,self.b,self.T,self.delta,self.M,self.g,self.f,self.finverse = u0,a,b,T,delta,M,g,f,finverse


		self.u = [[]]
		self.pos = [[]]

		N = int(self.M*(self.b-self.a)/self.delta)
		dx = (self.b-self.a)/N
		for j in range(0,N-1,1):
			self.u[0].append((1/dx)*integrate.quad(self.u0,self.a+j*dx,self.a+(j+1)*dx)[0]) 
			self.pos[0].append(self.a+(j+1)*dx)
		self.u[0].append((1/dx)*integrate.quad(self.u0,self.a+(N-1)*dx,self.b)[0])

		self.s = [[]]
		self.collisions = [0.]

		self.absolute_collision_times = [0.]


		# self.u = [[2.,-1.,1.]]
		# self.pos = [[-2.,0.]]
		# print self.IteratedRiemann(self.u[0],self.pos[0])
		# print self.CollisionTime(self.IteratedRiemann(self.u[0],self.pos[0])[0],self.IteratedRiemann(self.u[0],self.pos[0])[1])

	def solveFrontTracking(self):
		T, delta, u, pos, s, collisions = self.T, self.delta, self.u, self.pos, self.s, self.collisions

		time_elapsed = 0.
		
		while time_elapsed <T:
			pos.append([])
			u.append([])
			### Advance positions to t=time_elapsed
			pos[-1], u[-1] = self.AdvancePositions(collisions[-1],pos[-2],s[-1],u[-2])
			### Solve Riemann problem at t=time_elapsed
			if s[-1] != []:
				s.append([])
			s[-1], pos[-1], u[-1] = self.IteratedRiemann(u[-1],pos[-1])

			collisions.append(self.CollisionTime(s[-1],pos[-1],u[-1]))
			if collisions[-1] == 0.:
				break
			elif collisions[-1] < 0.:
				print 'Error: Negative collision time'
				break
			else:
				time_elapsed += collisions[-1]
		pos.append([])
		u.append([])
		pos[-1], u[-1] = self.AdvancePositions(T-time_elapsed+collisions[-1],pos[-2],s[-1],u[-2])

		for i in range(0,len(collisions)-2,1):
			self.absolute_collision_times.append(self.absolute_collision_times[i] + collisions[i+1])
		self.absolute_collision_times.append(self.T)

		# print '_______________________________________________________________________________________________'
		# print '_______________________________________________________________________________________________'
		# print 'coll = %s' %collisions
		# print '_______________________________________________________________________________________________'
		# print 'abs_coll = %s' %self.absolute_collision_times
		# print '_______________________________________________________________________________________________'
		# print 's = %s' %s
		# print '_______________________________________________________________________________________________'
		# print 'pos = %s' %pos
		# print '_______________________________________________________________________________________________'
		# print 'u = %s' %u
		# print '_______________________________________________________________________________________________'
		# print '_______________________________________________________________________________________________'

	def IteratedRiemann(self,u,pos):
		N = len(u)
		if N != len(pos)+1:
			print 'wrong dimensions'
		s = []
		pos_new = []
		u_new = [self.u[0][0]]
		for i in range(0,N-1,1):
			s.extend(self.SingleRiemann(u[i],u[i+1],pos[i])[0])
			pos_new.extend(self.SingleRiemann(u[i],u[i+1],pos[i])[1])
			u_new.extend(self.SingleRiemann(u[i],u[i+1],pos[i])[2])
		return s, pos_new, u_new

	def checkEqual(self,iterator):
   		return len(set(iterator)) <= 1

	def CollisionTime(self,s,pos,u):
		### insert x=0 into pos and s=0 into s
		if self.f != self.g:
			# print 'Yeah, we have a flux discontinuity'
			for i in range(0,len(s),1):
				if pos[i] == 0.:
					break
				if pos[i] > 0.:
					print 'Insert discontinuity at 0'
					print pos[i]
					pos.insert(i,0.)
					s.insert(i,0.)
					u.insert(i+1,u[i])
					break
				if i == len(s)-1:
					pos.append(0.)
					s.append(0.)
					u.append(u[-1])

		### calculate time until next collision
		if len(s) <=1 or self.checkEqual(s):
			return self.T + 1
		# print pos
		t_coll = []
		for i in range(0,len(s)-1,1):
			s_curr = s[i]
			pos_curr = pos[i]
			s_next = s[i+1]
			pos_next = pos[i+1]

			if s_curr <= s_next:
				continue
			elif pos_curr == pos_next and s_curr != s_next:
				# t_coll.append(0.)
				t_coll.append(self.T + 1)
			else:
				t_coll.append((pos_next - pos_curr)/(s_curr - s_next))
				# if (pos_next - pos_curr)/(s_curr - s_next) < 0:
					# print 'We have a negative collision time:'
					# print pos_next - pos_curr
					# print 'pos_curr %.16f' %pos_curr
					# print 'pos_next %.16f' %pos_next
					# print 's_curr %f' %s_curr
					# print 's_next %f' %s_next
		if t_coll == []:
			return self.T + 1.
		# print pos
		# print t_coll
		return min(t_coll)

	def AdvancePositions(self,collisionTime,pos,s,u):
		if collisionTime == 0.:
			return pos, u
		pos_new = []
		u_new = list(u)
		# print pos

		for i in range(0,len(pos),1):
			pos_new.append(pos[i]+s[i]*collisionTime)

		i = 0
		while i < len(pos_new)-1:
			# if pos_new[i] == pos_new[i+1]:
			if np.isclose(pos_new[i],pos_new[i+1]):
				pos_new.pop(i)
				u_new.pop(i+1)
				i-=1
			i+=1
		
		return pos_new, u_new

	def SingleRiemann(self,uL,uR,xpos):
		u_new = []
		s_new = []
		pos_new = []
		uM = uL
		delta = self.delta

		if xpos == 0.:
		# if np.isclose(xpos,0):
			uM = self.finverse(self.g(uL))
			############## Remove this part #########################################
			if uM == -uL:
				uM *=-1.
			############## Remove this part #########################################
			if uM != uL:
			# if not np.isclose(uM,uL):
				u_new.append(uM)
				s_new.append(0.)
				pos_new.append(0.)

		if xpos < 0:
			flux = self.g
		if xpos >= 0:
			flux = self.f
		### solve standard Riemann problem between uM and uR for convex flux
		if uM > uR:
			j = int(np.floor(uM/delta))
			k = int(np.floor(uR/delta))
			# s_new.append((flux(uM)-flux(uR)/(uM-uR))) ### This is wrong
			s_new.append((flux(k*delta)+(uR-k*delta)*(flux((k+1)*delta)-flux(k*delta))/delta - flux(j*delta) - (uM-j*delta)*(flux((j+1)*delta)-flux(j*delta))/delta)/(uR-uM))
			u_new.append(uR)
			pos_new.append(xpos)
			return s_new, pos_new , u_new
		elif uM < uR:
			### get index j such that uM in [(j-1)delta, j delta)
			j = int(np.ceil(uM/delta))
			if j*delta == uM:
				j += 1
			### get index k such that uR in [k delta, (k+1)delta)
			k = int(np.floor(uR/delta))
			if k*delta == uR:
				k -= 1

			s_new.append((flux(j*delta) - flux((j-1)*delta))/delta)
			pos_new.append(xpos)
			for i in range(j,k+1,1):
				s_new.append((flux((i+1)*delta)-flux(i*delta))/delta)
				pos_new.append(xpos)
				u_new.append(i*delta)
			u_new.append(uR)
			return s_new, pos_new, u_new
		elif uM == uR and xpos == 0.:
			return [0.], [0.], [uR]
		return [], [], []








""" Some example fluxes and initial data """

def Burgers(u):
	return 0.5*u**2
def BurgersInverse(u):
	return np.sqrt(2*u)
def Transport(u):
	return u
def BuckleyLeverett(u):
	return (u**2)/(u**2+(1-u)**2)
def BuckleyLeverettInverse(u):
	if u == 0.5:
		return 0.5
	else:
		return (u-np.sqrt((1-u)*u))/(2*u-1)
def SlowBuckleyLeverett(u):
	return 0.8*(u**2)/(u**2+(1-u)**2)





def u0test(x):
	if x < -2.:
		return 2.
	elif x < 0:
		return -1.
	else:
		return 1.
def u0Riemann(x):
	if x < -0.5:
		return 0.5
	else:
		return 2.
def u0bump(x):
	return 2 + np.exp(-100*(x+0.75)**2)
def u0downstep(x):
	if x < 0:
		return 1.
	else:
		return 0.
def u0podium(x):
	### Note that the code only works for convex flux atm
	if (x>=-1) and (x<=0):
		return 0.4
	else:
		return 0.2



def Experiment(u0,delta,M):
	""" Takes the initial datum, delta and M and sets up the corresponding experiment """
	if u0 == u0Riemann:
		### Switch from Transport to Burgers across x=0
		a = -1.
		b = 1.
		T = 0.9
		g = Transport
		f = Burgers
		finverse = BurgersInverse
	elif u0 == u0bump:
		### Switch from Burgers to Transport across x=0
		a = -1.
		b = 1.
		T = 0.5
		g = Burgers
		f = Transport
		finverse = Transport
	elif u0 == u0podium:
		a = -1.1
		b = 1.6
		T = 1.
		g = SlowBuckleyLeverett
		f = BuckleyLeverett
		finverse = BuckleyLeverettInverse
	else:
		### Burgers flux on both sides of x=0
		a = -3.
		b = 6.
		T = 7.
		g = Burgers
		f = Burgers
		finverse = BurgersInverse
		
	method = Method(u0=u0, a=a, b=b, T=T, delta=delta, M=M, g=g, f=f, finverse=finverse)
	method.solveFrontTracking()
	return method






pause = True

def plot(method):
	""" Plots an animation of the solution in time generated by Front Tracking """

	a, b, u0, T, delta, M = method.a, method.b, method.u0, method.T, method.delta, method.M
	u0vec = np.vectorize(u0)
	xRef = np.linspace(a,b,512)
	
	if u0 == u0bump:
		def helpfunctional(eta,x):
			return x-u0(eta)*T - eta
		def Uexact(x):
			eta0 = brentq(helpfunctional,-100.,100.,args=(x))
			if eta0 < a:
				return u0(a)
			elif eta0 > b:
				return u0(b)
			else:
				return u0(eta0)
		uexactvec = np.vectorize(Uexact)
	
	x1 = [a]
	x1.extend(method.pos[-1])
	y1 = list(method.u[-1])
	
	if method.pos[-1][-1] < b:
		x1.append(b)
		y1.append(method.u[-1][-1])
	
	
	
	def onClick(event):
		global pause
		pause ^= True
	def onClose(event):
		global pause
		pause = True
	
	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	
	if u0 == u0Riemann:
		ymin = 0.
		ymax = 2.5
		t_increment = 0.01
	elif u0 == u0bump:
		ymin = 1.
		ymax = 4.
		t_increment = 0.01
	elif u0 == u0podium:
		ymin = 0.
		ymax = 1.
		t_increment = 0.01
	else:
		ymin = -2.
		ymax = 4.
		t_increment = 0.05
	ax1.axis([a,b,ymin,ymax])
	
	line1, = ax1.step([],[], where='post')
	line2, = ax1.plot([],[],'r')
	line3, = ax1.plot([],[],'g--')
	line4, = ax1.step([],[],'r--', where='post')
	
	lines = [line1, line2, line3, line4,]
	
	def simData():
		t = 0.
	
		x = [a]
		x.extend(method.pos[0])
		x.append(b)
	
		y = list(method.u[0])
		y.append(method.u[0][-1])
		
		absolute_time_of_next_collision = 0.
		# yield x, y
		for i in range(0,len(method.collisions)-1,1):
			absolute_time_of_last_collision = absolute_time_of_next_collision
			if method.collisions[i+1] == 0.:
				absolute_time_of_next_collision = T
			else:
				absolute_time_of_next_collision += method.collisions[i+1]
			while t < min(absolute_time_of_next_collision, T):
				if not pause:
					x = [a]
					for j in range(0,len(method.pos[i+1]),1):
						x.append(method.pos[i+1][j] + (t-absolute_time_of_last_collision)*method.s[i][j])
	
					y = list(method.u[i+1])
	
					if method.pos[i+1][-1] < b:
						x.append(b)
						y.append(method.u[i+1][-1])
	
					t += t_increment
				yield x, y
	
	def frame(simData):
		x, y = simData[0], simData[1]
		lines[0].set_data(x,y)
	
		# if u0 == u0bump:
		# 	lines[1].set_data(xRef,uexactvec(xRef))
		lines[2].set_data(xRef, u0vec(xRef))
		lines[3].set_data(x1,y1)
		return lines
	
	fig.canvas.mpl_connect('button_press_event', onClick)
	fig.canvas.mpl_connect('close_event', onClose)
	
	anim = animation.FuncAnimation(fig, frame, simData, interval = 20, blit=False,repeat=False)
	
	major_ticks = np.arange(ymin, ymax, delta)
	ticks2 = np.arange(a, b, delta/M)
	ax1.set_yticks(major_ticks,minor=True)
	ax1.set_xticks(ticks2,minor=True)
	ax1.grid(which='minor', alpha=1.)
	plt.show()




def L1Error(method1, method2):
	""" Calculates the L1 error between two Front Tracking solutions """
	a, b = method1.a, method1.b
	if (a != method2.a) or (b!= method2.b):
		print 'Wrong interval lengths'
	pos1 = list(method1.pos[-1])
	pos1.append(b)

	pos2 = list(method2.pos[-1])
	pos2.append(b)

	u1, u2 = list(method1.u[-1]), list(method2.u[-1])



	### delete the part of the solution that is not in [a,b]
	while pos1[-2] >= pos1[-1]:
		pos1.pop(-2)
		u1.pop(-1)
	while pos2[-2] >= pos1[-1]:
		pos2.pop(-2)
		u2.pop(-1)

	### check whether the positions are sorted
	if pos1 != sorted(pos1):
		print 'pos1 not sorted'
	if pos2 != sorted(pos2):
		print 'pos2 not sorted'

	### sanity check
	if not ((len(pos1) == len(u1)) and (len(pos2) == len(u2))):
		print 'Error in solution length'

	i, j = 0, 0
	error = 0.
	x_prev = a
	x_curr = a
	while (i < len(u1)) and (j < len(u2)):
		### get x_curr
		x_curr = min(pos1[0],pos2[0])

		error += (x_curr - x_prev)*np.abs(u1[i] - u2[j])

		x_prev = x_curr
		### increment i or j
		if x_curr == pos1[0]:
			i += 1
			pos1.pop(0)
		if x_curr == pos2[0]:
			j += 1
			pos2.pop(0)
	return error



def ConvergenceRates(u0):
	""" 
	Saves a table containing the convergence rates of the experiment corresponding to u0 in {u0Riemann,u0bump}
	as a txt file 
	"""
	Nmin = 16
	Nmax = 1024
	NFine = 2048

	methodFine = Experiment(u0=u0, delta=1./NFine, M=1.)

	tablename = 'Table'
	tablename += str(NFine)
	if u0 == u0Riemann:
		tablename += 'RiemannID_Transport_to_Burgers'
	elif u0 == u0bump:
		tablename += 'BumpID_Burgers_to_Transport'
	tablename += '.txt'

	table = ''
	k = int(np.log2(Nmax/Nmin))
	err = np.zeros(k+1)
	rate = np.zeros (k+1)

	for i in range(0,k+1,1):
		Ncomp = Nmin*(2**i)
		methodcomp = Experiment(u0=u0, delta=1./Ncomp, M=1.)

		err[i] = L1Error(methodcomp,methodFine)
		if i>0:
			rate[i] = -np.log2(err[i]/err[i-1])
		table += "N = %i L1 error: %.3e Rate %.2f \n" %(Ncomp,err[i],rate[i])
	np.savetxt(tablename, [table],fmt='%s')


def saveSolution(method):
	text = 'Solution'
	if method.u0 == u0Riemann:
		text += 'RiemannID_Transport_to_Burgers'
	elif method.u0 == u0bump:
		text += 'BumpID_Burgers_to_Transport'
	text += str(method.T)
	text += '.txt'

	x = list(method.pos[-1])
	x.insert(0,method.a)
	help = ''
	for j in range(0,len(method.u[-1]),1):
		help += "%.4f %.4f \n" %(x[j], method.u[-1][j])

	help += "%.4f %.4f \n" %(method.b, method.u[-1][-1])
	help += "\n"

	np.savetxt(text, [help], fmt='%s')

def saveExperiments(u0, delta, M):
	t = [0., 0., 0.]
	if u0 == u0Riemann:
		a = -1.
		b = 1.
		g = Transport
		f = Burgers
		finverse = BurgersInverse
		t[0] = 0.3
		t[1] = 0.6
		t[2] = 0.9
	elif u0 == u0bump:
		a = -1.
		b = 1.
		g = Burgers
		f = Transport
		finverse = Transport
		t[0] = 0.2
		t[1] = 0.3
		t[2] = 0.5
	for i in range(0,3,1):
		method = Method(u0=u0, a=a, b=b, T=t[i], delta=delta, M=M, g=g, f=f, finverse=finverse)
		method.solveFrontTracking()
		saveSolution(method)


def saveID(u0):
	help = ''
	if u0 == u0Riemann:
		text = 'RiemannID.txt'
		help += "-1. 0.5 \n"
		help += "-0.5 2. \n"
		help += "1. 2. \n"
		help += "\n"
	elif u0 == u0bump:
		text = 'BumpID.txt'
		u0vec = np.vectorize(u0)
		x = np.linspace(-1.,1.,512)
		y = u0vec(x)
		for i in range(0,len(x),1):
			help += "%.4f %.4f \n" %(x[i],y[i])
		help += "\n"


	np.savetxt(text, [help], fmt='%s')



if __name__ == '__main__':
	### Experiment 1
	method1 = Experiment(u0=u0Riemann, delta=1./32, M=1.)

	### Experiment 2
	method2 = Experiment(u0=u0bump, delta=1./32, M=1.)

	plot(method1)
	plot(method2)

	# method3 = Experiment(u0=u0podium, delta=1./10, M=1.)
	# plot(method3)


	### Convergence Rates for Experiment 1
	# ConvergenceRates(u0Riemann)

	### Convergence Rates for Experiment 2 (only uncomment if you have some time on your hand)
	# ConvergenceRates(u0bump)








