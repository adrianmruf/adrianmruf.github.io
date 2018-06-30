import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import lines



class Problem(object):
	def __init__(self, alpha, beta, gamma, T):
		self.alpha, self.beta, self.gamma, self.T = alpha,beta,gamma,T



class Method(object):
	def __init__(self,problem,u0,M,N,bc,fluxind):
		self.problem=problem
		self.M, self.N = M, N
		self.u0 = u0
		self.bc = bc
		self.fluxind = fluxind

		self.u=np.zeros((M+2,N+1))

		N = self.N
		dx = 1./N

		self.u[0,0] = 2*N*integrate(0,0.5*dx,u0,5)
		self.u[0,N] = 2*N*integrate((N-0.5)*dx,1,u0,5)
		for j in range(1,N,1):
			self.u[0,j] = N*integrate((j-0.5)*dx,(j+0.5)*dx,u0,10)
			# self.u[0,j]=self.u0(j*dx)

	def solve(self):
		N = self.N
		M = self.M
		T = self.problem.T
		dx = 1./N
		dt = float(T)/(M+1)
		lam = dt/dx
		gamma = self.problem.gamma

		u = self.u


		for n in range(0,M+1,1):
			F = self.Flux(u[n,:])
			P = self.P(u[n,:])

			u[n+1,:] = u[n,:] - lam*(F[1:] - F[:-1]) + gamma*dt*P
			
			if (self.fluxind != 'Lax-Friedrichs') and (self.fluxind != 'Engquist-Osher'): 
				print 'wrong flux indication'
				break
			if self.bc == 'dirichlet':
				u[n+1,0] = 1/dt*integrate(n*dt,(n+1)*dt,self.problem.alpha,10)
				u[n+1,N] = 1/dt*integrate(n*dt,(n+1)*dt,self.problem.beta,10)
				# u[n+1,0] = self.problem.alpha((n+1)*dt)
				# u[n+1,N] = self.problem.beta((n+1)*dt)
			if (self.bc != 'dirichlet') and (self.bc != 'periodic'):
				print 'wrong bc indication'
				break

			


	def Flux(self,u):
		N = self.N
		M = self.M
		T = self.problem.T
		dx = 1./N
		dt = float(T)/(M+1)

		F = np.zeros(N+2)
		for j in range(0,N,1):
			ul = u[j]
			ur = u[j+1]
			
			if self.fluxind == 'Lax-Friedrichs':
				# uncomment if the linear flux f(u)=u is desired
				# F[j+1] = 0.5*((ul+ur) + dx/dt*(ul-ur))

				# uncomment if the flux f(u)=0.5*u**2 is desired
				F[j+1] = 0.5*(0.5*(ul**2+ur**2) + dx/dt*(ul-ur))

			if self.fluxind == 'Engquist-Osher':
				if ul>0 and ur>=0:
					F[j+1] = .5*ul**2
				if ul<=0 and ur<0:
					F[j+1] = .5*ur**2
				if ul<0 and ur>0:
					F[j+1]=0
				if ul>0 and ur<0:
					F[j+1] = .5*(ul**2+ur**2)


		# F has the form [0, F(u0,u1), ... , F(uN,uN-1),0] now
		# this is only necessary for periodic boundary conditions but should not make a difference w/o bc
		F[0] = F[N]
		F[N+1] = F[1]
		return F

	def P(self,u):
		N = self.N
		dx = 1./N
		P1 = [0.]
		P1.extend(dx*np.cumsum(u[:-1]))
		P1 = np.asarray(P1)
		P = P1 + 0.5*dx*u

		if self.bc == 'dirichlet':
		# 	return P
		# if self.bc == 'periodic':
			return P - dx*np.sum(P[:-1])


def U0trigon(x):
	return -0.05*np.cos(2*np.pi*x)

def U0cornerwave(x):
	# return 0
	if 0<=x<=0.5:
		return 1./6*(x-0.5)**2 + 1./6*(x-0.5) + 1./36
	elif 0.5<x<=1:
		return 1./6*(x-0.5)**2 - 1./6*(x-0.5) + 1./36
	else:
		print 'error'
		return -1/36.

# Integration routine
def integrate(a,b,f,N):
    s = 0.
    dx = float(b-a)/N
    for i in range(N):
        s += f(a+i*dx)
    return s * dx


# This is the explicit solution for the cornerwave
def Uexp(x,t):
	if 0<=x-t/36.<=1:
		return U0cornerwave(x-t/36.)
	elif -1<=x-t/36.<0:
		return U0cornerwave(1+x-t/36.)
# Dirichlet bc can only be used for the cornerwave because there the explicit solution is known

def alpha(t):
	# return 0
	return Uexp(0,t)
def beta(t):
	return 0
	return Uexp(1,t)






problem = Problem(alpha,beta,gamma=1,T=36)

###
# N=64
# M=91
###
N=128
M=183
###
# N=256
# M=367
###
# N=512
# M=736
###
# N=1024
# M=1473
###
# flux = 'Lax-Friedrichs'
flux = 'Engquist-Osher'

method = Method(problem, u0=U0cornerwave, M=M, N=N, bc='dirichlet', fluxind=flux)
method.solve()




#####
# Calculate the L1 error between approximation and exact solution (Cornerwave)
#####

# err = 0
# dx = 1./method.N

# def errfct(x):
# 		return np.abs(method.u[method.M+1,0]-U0cornerwave(x))

# err += integrate(0,0.5*dx,errfct,5)

# def errfct(x):
# 		return np.abs(method.u[method.M+1,method.N]-U0cornerwave(x))

# err += integrate((method.N-0.5)*dx,1,errfct,5)

# for j in range(1,method.N,1):

# 	def errfct(x):
# 		return np.abs(method.u[method.M+1,j]-U0cornerwave(x))

# 	err += integrate((j-0.5)*dx,(j+0.5)*dx,errfct,10)

# print "L1 error:" ,err


#######
# Calculate error between approximate solution and finer approximation
#######

Nmax = 2048
Mmax = 2948
method2 = Method(problem, u0=U0cornerwave, M=Mmax, N=Nmax, bc='dirichlet', fluxind='Engquist-Osher')
method2.solve()

# dx = 1./method.N
# dx2 = 1./method2.N
# ratio = method2.N/method.N

# u2new = np.zeros(2*method2.N)
# for j in range(1,method2.N,1):
# 	u2new[2*j-1] = method2.u[method2.M+1,j]
# 	u2new[2*j] = method2.u[method2.M+1,j]
# u2new[0] = method2.u[method2.M+1,0]
# u2new[-1] = method2.u[method2.M+1,-1]

# unew = np.zeros(2*method2.N)
# for j in range(0,ratio,1):
# 	unew[j] = method.u[method.M+1,0]
# 	unew[-j-1] = method.u[method.M+1,-1]
# for j in range(0,method.N-1,1):
# 	for i in range(0,2*ratio,1):
# 		unew[ratio+j*2*ratio+i] = method.u[method.M+1,j+1]


# err = 0
# for j in range(0,2*method2.N,1):
# 	err += 0.5*dx2*np.abs(u2new[j]-unew[j])

# print "L1 error between solutions on coarse and fine grid:" ,err







#######################
#
# Plotting the solution
#
#######################
x = np.linspace(0,1,method.N+1)
x2= np.linspace(0,1,method2.N+1)

pause = True

def onClick(event):
	global pause
	pause ^= True

fig = plt.figure()
ax1 = fig.add_subplot(111)
# lines = ax1.plot([],[])

#
line1, = ax1.plot([],[],'k')
line2, = ax1.plot([],[],'k')
line3, = ax1.plot([],[],'k')
line4, = ax1.plot([],[],'k')
line5, = ax1.plot([],[],'k--')
lines = [line1, line2, line3, line4, line5]
#

if method.u0 == U0cornerwave:
	plt.axis([x[0],x[-1],-0.02,0.03])
if method.u0 == U0trigon:
	plt.axis([x[0],x[-1],-0.06,0.08])


def simData():
	t_max = method.M+1
	t = 0
	y = method.u[0,:]
	while t<t_max:
		if not pause:
			y = method.u[t,:]
			t += 1
		yield y,t

def frame(simData):
	y,t = simData[0],simData[1]
	# if I want the integral
	# time_text.set_text(time_template%(np.trapz(y,x)))
	# if I want the time t
	time_text.set_text(time_template%(t*method.problem.T/(method.M+1),method.M,method.N))
	# lines[0].set_data(x2,method2.u[0,:])
	lines[0].set_data(x,method.u[t,:])
	# lines[4].set_data(x2,method2.u[method2.M+1,:])

	# lines[1].set_data(x,method.u[75,:])
	# lines[2].set_data(x,method.u[150,:])
	# lines[3].set_data(x,method.u[225,:])
	# lines[4].set_data(x,method.u[method.M+1,:])
	# lines[3].set_data(x2,method2.u[method2.M+1,:])
	return lines

time_template = 't = %.0f, M = %.0f, N = %.0f'
time_text = ax1.text(0.45, 1.04, '', transform=ax1.transAxes)

fig.canvas.mpl_connect('button_press_event',onClick)

anim = animation.FuncAnimation(fig, frame, simData, interval = 50, blit=False,repeat=False)

plt.show()

Writer = animation.writers['ffmpeg']
writer = Writer(fps=20,bitrate=-1)
# anim.save('periodic.mp4', writer=writer)





