import numpy
import pymbar
import matplotlib.pyplot as plt 

R = 1.9872036e-3                          # kcal/mol.K
atm2Pa = 101325                           # 1 atm = 101325 Pa
J2kcal = 1/4184
NA = 6.022140857e23                       # Avogadro number
m3_2_A3 = 1e30
P = 1*atm2Pa*J2kcal*NA/m3_2_A3      # pressure in units of energy/volume [kcal/mol.A³]


mvv2e = 2390.057364     # Da*A²/fs² to kcal/mol
kB = 8.31451E-7         # in Da*A²/(fs²*K)
Pconv = 1.6388244954E+8 # Da/(A*fs²) to atm


U_index = 5
V_index = 4
prop_index = 5    # sera usado ET para o calculo de CP a partir de H = ET + PV

#  0-Step    |   1-Temp    |  2-Press   |  3-v_dens |  4-Volume  |  5-PotEng | 6-KinEng | 7-TotEng | 8-E_vdwl | 9-E_coul | 10-E_pair | 11-bond 
# 12-E_angle |  13-E_dihed | 14-E_impro | 15-E_mol  | 16-E_long  | 17-E_tail
# Lembrando que: Ucf = Epot - Emol
#                Hcf = Epot - Emol + P_nom*V 

#===================================================================================================
#  Definition of auxiliary functions
#===================================================================================================
def logsumexp(x):
	xmax = numpy.amax(x)
	return xmax + numpy.log(sum(numpy.exp(x - xmax)))

def average(a,x):
	xmax = numpy.amax(x)
	return sum(a*numpy.exp(x - xmax))/sum(numpy.exp(x - xmax))

#===================================================================================================
#  Temperatures Reading
#=================================================================================================== 
filename = 'temp.inp'
print ("Reading %s..." % filename)
with open(filename,'r') as infile:
	elements = infile.readline().split()
	K = len(elements)
	temp = []
	for k in range(K):
		temp.append( [elements[k]] )
temp = numpy.array(temp, float)
beta = 1.0/(R*temp)
K  = len(beta)
infile.close()

#===================================================================================================
#  Manage log files
#===================================================================================================
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

nsnaps = 0
for ind in range(K):
	filename = str(ind+1) + '.log'
	nsnaps = max(nsnaps,file_len(filename))
nsnaps = int(nsnaps)

#===================================================================================================
#  Load equilibrated data
#===================================================================================================

u_kln = numpy.zeros([K,K,nsnaps], numpy.float64)
A_kn = numpy.zeros([K,nsnaps], numpy.float64)

n = numpy.zeros(K, int)
for ind in range(K):
	filename = str(ind+1) + '.log'
	print ("Reading metadata from %s..." % filename)
	infile = open(filename, 'r')
	lines = infile.readlines()
	for (i,line) in enumerate(lines):
		elements = line.split()
		H = float(elements[U_index]) + P*float(elements[V_index])
		for k in range(K):
			u_kln[ind,k,n[ind]] = beta[k]*H
		A_kn[ind,n[ind]] = float(elements[prop_index])*float(elements[prop_index]) + P*float(elements[V_index])*float(elements[prop_index])
		#A_kn[ind,n[ind]] = float(elements[prop_index])
		n[ind] += 1

for l in range(K):
	plt.hist(u_kln[l,l,0:n[l]],bins= 100,density = True, edgecolor = "none", label = str(l+1))
plt.legend()
plt.show()

print ("Running MBAR...")
MBAR = pymbar.mbar.MBAR(u_kln, n, verbose = True, relative_tolerance = 1.0e-10,
                        solver_protocol='default', initialize = 'BAR')

# Get matrix of dimensionless free energy differences
results = MBAR.compute_free_energy_differences(uncertainty_method='svd-ew',
                                                                  return_theta = True) 

f = results['Delta_f']
df = results['dDelta_f']
#print ("Free energies:")
#print (f)
print ("Uncertainties:")
print (df)

#print ("Computing expectations...")
compute1 = MBAR.compute_expectations(A_kn) #EA_k e dEA_k
compute2 = MBAR.compute_expectations(A_kn*A_kn) #EA2_k e dEA2_k

EA_k = compute1['mu']
dEA_k =compute1['sigma']
EA2_k = compute2['mu']
dEA2_k =compute2['sigma']

print ("Expectations:")
print (EA_k)

print("Uncertainties:")
print (dEA_k)

ntotal = sum(n)
u0 = numpy.zeros(ntotal, float)
PE = numpy.zeros(ntotal, float)
A = numpy.zeros(ntotal, float)
i = 0
for j in range(K):
	for k in range(n[j]):
		u0[i] = -logsumexp(-u_kln[j,:,k] + f[0] + numpy.log(n))
		PE[i] = u_kln[j,j,k]/beta[j]
		A[i] = A_kn[j,k]
		i += 1
A2 = A*A

npoints = 1000   # numero de pontos que serão formados a partir do tratamento dos dados
Tmin = min(temp)
Tmax = max(temp)
T_curve = numpy.zeros(npoints, float)
f_curve = numpy.zeros(npoints, float)
A_curve = numpy.zeros(npoints, float)
A2_curve = numpy.zeros(npoints, float)
for point in range(npoints):
	T_curve[point] = point*(Tmax - Tmin)/npoints + Tmin
	du = u0 - PE/(R*T_curve[point])
	f_curve[point] = -logsumexp( du )
	A_curve[point] = average(A, du + f_curve[point] )
	A2_curve[point] = average(A2, du + f_curve[point] )

Var = (EA2_k - EA_k*EA_k)/EA_k
dVar = numpy.zeros(K, float)

#fig, ax = plt.subplots(3,1,sharex=True)
fig, ax = plt.subplots(2,1,sharex=True)

#ax[0].errorbar(temp, f, yerr=df, fmt='o')
#ax[0].plot( T_curve, f_curve )

ax[0].errorbar(temp, EA_k, yerr=dEA_k, fmt='o')
ax[0].plot(T_curve, A_curve )

ax[1].errorbar(temp, Var, yerr=dVar, fmt='o')
ax[1].plot(T_curve, (A2_curve - A_curve*A_curve)/A_curve,'.')
plt.show()

numpy.savetxt('resultUH.txt', A_curve)
numpy.savetxt('result2.txt', Var)
numpy.savetxt('temp.txt', temp)
numpy.savetxt('T.txt', T_curve)
