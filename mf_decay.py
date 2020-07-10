import sys
from math import *
import numpy as np
import matplotlib as mpl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Palatino']})
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
mpl.rcParams.update({'font.size': 15})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 14})

def t_cascade (B0):

	res =  3.0 * 1e4 * 1e15 / B0
	return res


def temp (t, B0):
	T = 7e8 *exp(-t/150.0) + 1.5e8 * exp(-t/2.5e6) ## Additional heating of magnetars which cannot occur due to normal NS cooling

	if t > t_cascade (B0):

		T = 6.56e8 * pow (t, -0.185) * exp(-t / 8.58e5)


	return T

#def temp_surf (t):
#
#        T = 2.1938e6 * pow (t, -0.0964) * exp(-t/1.5285e6)
#        return T


def tau_ohm_ph (t, B0):

	T_curr = temp (t, B0)

#	print ('Age: t = ', t, ' T_curr = ', T_curr)

	if (T_curr > 0):
		res = 2e6 * pow(T_curr/1e8, -2.0)
	else:
		res = 1e12

	return res

def tau_ohm_q (t, Q):

	res = 2e6 / Q

	return res


def tau_ohm(t, Q, B0):

	t1_curr = tau_ohm_ph (t, B0)
	t2_curr = tau_ohm_q (t, Q)

	res = 1.0 / t1_curr  + 1.0 / t2_curr

	res = 1.0 / res

	return res


def tau_hall (t, B0):

	t_casc = t_cascade (B0)	

	if (t < t_casc):
		res = 1e4 * 1e15 / B0 
	else:
		res = 1e12

	return res

def tau_hall_adapted (t, B, B0):

	t_casc = t_cascade (B0)	

	if (t < t_casc):
		res = 1e4 * 1e15 / B 
	else:
		res = 1e12

	return res


def deq_rhs (t, Q, B0, B):

	tau_ohm_  = tau_ohm  (t, Q, B0)
	tau_hall_ = tau_hall (t, B0) 

	res = - B / tau_ohm_ - B*B / (B0 * tau_hall_)
	return res

def actual_B1 (t, Q, B0):

	delta_t = 100.0

	n = int (t / delta_t)

	B = B0

	list_B = []

	for i in range (0, n):

		list_B.append(B)

		t_curr = i * delta_t

		k1 = delta_t * deq_rhs (t_curr, Q, B0, B)
		k2 = delta_t * deq_rhs (t_curr + 0.5 * delta_t, Q, B0, B + 0.5 * k1)
		k3 = delta_t * deq_rhs (t_curr + 0.5 * delta_t, Q, B0, B + 0.5 * k2)
		k4 = delta_t * deq_rhs (t_curr + delta_t, Q, B0, B + k3)
		
		B_new = B + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
		B = B_new

	return list_B
def get_P (t, B_list, P0, Q):

	delta_t = 100.0

	beta = 1.6e-39 ## G^2 / s

	n = int (t / delta_t)

	P_list      = []
	dotP_list   = []
	brkInd_list = []
	tau_list    = []

	P_list.append(P0)
	dotP_new = 2.0 / 3.0 * beta * pow(B_list[0], 2.0) / P0
	dotP_list.append(dotP_new)
	#brkInd_list.append(3.0)

	for i in range (0, n):

		t_curr = i * delta_t		

		P_new      = sqrt(pow(P_list[-1], 2.0) + 2.0 / 3.0 * delta_t * 3.154e+7 * pow(B_list[i], 2.0) * beta)

		dotP_new   = 2.0 / 3.0 * beta * pow(B_list[i], 2.0) / P_new

		brkInd_new = 3.0 - 2.0 * P_new / dotP_new * deq_rhs (t_curr, Q, B_list[0], B_list[i]) / B_list[i] / 3.154e+7

		tau_new = P_new / (2.0 * dotP_new) / 3.154e+7

		P_list.append(P_new)
	
		dotP_list.append(dotP_new)

		brkInd_list.append(brkInd_new)

		tau_list.append(tau_new)

	return [P_list, dotP_list, brkInd_list, tau_list]

t = 1e7
Q = 10
B0 = 1e15
P0 = 0.04

val1 =  actual_B1  (t, Q, B0)

val_t = []
val_time = []

for i in range (0, 30):

	val_time.append (t/30*i)
	val_t.append (temp (t/30*i, B0))

	print (i, val1[i])
	
#plt.plot (val_time, val_t)
#plt.xlabel('Age (yrs)')
#plt.ylabel('T (K)')
#plt.xscale('log')
#plt.yscale('log')
#plt.show()


t_val = np.linspace (0, t, int(t/100.0))


plt.plot (t_val, val1, 'r-',  label='Magnetic field evolution')

plt.xlabel('T (yrs)')
plt.ylabel('B (G)')

plt.xscale('log')
plt.yscale('log')
plt.ylim([1e12,2e15])
plt.legend()
plt.savefig('different_approach.pdf')
plt.tight_layout()
plt.savefig('B_p.pdf')
plt.show()


sys.exit(0)

P_list, dotP_list, brkInd_list, tau_list = get_P (t, val1, P0, Q)

for i in range (0, len(P_list) - 1):

	print (i, P_list[i], dotP_list[i], brkInd_list[i])


plt.scatter (P_list, dotP_list)
plt.xlim([0.01, 1])
plt.ylim([1e-18, 1e-8])
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$P$ (sec)')
plt.ylabel(r'$\dot P$')
plt.show()

plt.plot (t_val, brkInd_list)
plt.xscale('log')
plt.xlabel('Age (years)')
plt.ylabel(r'$n$')
plt.show()

plt.plot (tau_list, brkInd_list)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\tau$ (years)')
plt.ylabel(r'$n$')
plt.savefig('n_of_tau.pdf')
plt.show()


f = open('all_charact.txt', 'w')

f.write('Age \t B \t P \t dotP \t n \t tau \n')

for i in range (0, len(t_val)):

	f.write(str(t_val[i]) + "\t" + str(val1[i]) + '\t' + str(P_list[i]) + '\t' + str(dotP_list[i]) + '\t' + str(brkInd_list[i]) + '\t' + str(tau_list[i]) + '\n')





