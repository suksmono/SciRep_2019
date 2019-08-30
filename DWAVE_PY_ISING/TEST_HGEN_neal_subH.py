"""
-------------------------------------------------------------------------------
Created on Wed Jul 18 2018
author: Suksmono@{STEI-ITB, MDR Inc.}
Find a Hadamard matrix using QA Computer/DWave
> Sub-problem: orhtogonalization of a set of vector
> N: subset cardinality, M: order of H-matrix
-------------------------------------------------------------------------------
"""
from sympy import *
from prb_symbolic import *
import neal
import numpy as np

# define matrix order
M=12 # iorder of matrix, dimension of column vetor
N=3 # cardinality or the number members in ortho-set
# number of required qubits
NQ=M*N + M*int(N*(N-1)/2)
'''
-------------------------------------------------------------------------------
1. Formulate Hks
-------------------------------------------------------------------------------
'''
"""
for 4x4 case, the arrangement of the qubits is
-------------------------------------------------------------------------------
 |<---problem--->|<------ ancillas ------>|
-------------------------------------------------------------------------------
             s-DOMAIN
  s0 s4 s8  | s12 s16 s20 
  s1 s5 s9  | s13 s17 s21 
  s2 s6 s10 | s14 s18 s22 
  s3 s7 s11 | s15 s19 s23 
-------------------------------------------------------------------------------
"""
ss=symbols('s0:%d'%NQ)
ss0=np.asarray(ss)
s=ss0.reshape(int(NQ/M),M).tolist()

print('Calculating Hks ...')

'''
Hks=0;
for m in range(0,N):
    ts=s[m]
    for nn in range(m+1,N):
        Hks = Hks + (np.dot(s[m],s[nn]))**2 
Hks=expand(Hks)
'''

Hks=formHks(N,s)

'''
---------------------------------------------------
simplify by substition of all si**2 terms: si**2->1
---------------------------------------------------
'''
print('Substitution of si**2->1 ...')
Hks=rmvIdSquare(Hks,s)

'''
-------------------------------------------------------------------------------
2. Transform Hks -> Hkq
-------------------------------------------------------------------------------
'''

"""
-------------------------------------------------------------------------------
 |<---problem--->|<------ ancillas ------>|
-------------------------------------------------------------------------------
             q-DOMAIN
  q0 q4 q8  | q12  q16 q20 
  q1 q5 q9  | q13  q17 q21 
  q2 q6 q10 | q14  q18 q22  
  q3 q7 q11 | q15  q19 q23
-------------------------------------------------------------------------------
"""
qq=symbols('q0:%d'%NQ)
qq0=np.asarray(qq)
q=qq0.reshape(int(NQ/M),M).tolist()

print('Transform: Hks->Hkq ...')
Hkq=Hks2Hkq(Hks, s, q)

'''
-------------------------------------------------------------------------------
2. Transform Hkq -> H2q -> H2s
-------------------------------------------------------------------------------
'''

'''
---------------------------------------------------
 define rows of substitution pair [i,j,k]: qi*qj->qk
---------------------------------------------------
'''

spair=genSubsPair(N)

'''
spair= [ [0,1,4], [0,2,5], [0,3,6], \
         [1,2,7], [1,3,8], [2,3,9] ]
'''
#
[NPAIR, XX]=np.shape(spair)
dijMax=M*M #**2
#HMax=NPAIR*dijMax
HMax=dijMax
delta=5*HMax #6*16 #2*NPAIR*M #2*(M**2)
print('Transform: Hkq->H2q->H2s ...')
H2q, H2s = q2s_symbolig(Hkq, q, s, spair, delta)
print('H2q:\n', H2q)
print('H2s:\n', H2s)

'''
------------------------------------------------------------
3. EXTRACT ISING PARAMETERS FROM SYMBOLIC SOLUTION
------------------------------------------------------------
'''
print('Obtaining Ising coefficients ...')
b, hi, Jij = isingCoeffs(H2s,NQ)


# normalize coefficients
maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij))])
hi=hi/maxCoeff
Jij=Jij/maxCoeff
#
b=b/maxCoeff
'''
-----------------------------------------------------------------------------
convert the problem into Ising coefficients
-----------------------------------------------------------------------------
'''
#in dictionary format
h={0:0}
J={(0,1):1}

for m in range(0,len(hi)):
    h[m]=hi[m]
    for n in range (m+1,len(hi)):
        J[m,n]=Jij[m,n]
    
'''
-----------------------------------------------------------------------------
4. SOLVE THE PROBLEM
-----------------------------------------------------------------------------
select a solver
> dimod: ExaxtSolver
> neal:  SimulatedAnnealingSampler
'''
#
print('Solving the problem using neal  ...')
#### add timing
import time
t0 = time.time()
##
solver=neal.SimulatedAnnealingSampler()
NSWEEPS=1*5*10*10*1000
response = solver.sample_ising(h, J, sweeps=NSWEEPS, num_reads=10)
#
elapsed = time.time()-t0
print('(simulated) annealing time = ', elapsed)
#vE=response.data_vectors['energy']
#aSol=response.samples_matrix

vE=response.record['energy']
aSol=response.record['sample']

#
print('Configurations:\n', aSol)
print('Energy:\n',vE)
#
idxMinE=np.argmin(vE)
print('Minimum Energy:',vE[idxMinE], 'supposed to be', -b)
print('Minimum Configurations:',aSol[idxMinE])
tSol=aSol[idxMinE]
vSol=tSol #[0]
#
HNM=vSol[0:M*N].reshape(N,M)
print('Found matrix:\n',HNM)
print('Indikator matrix:\n', \
       np.matmul(HNM.tolist(), HNM.transpose().tolist() ))
print('Subset cardinality:',N, ', H-Order', M)
print('Number of qubits:', NQ,', number of H2s terms:',\
      len(H2s.as_coefficients_dict()) )
