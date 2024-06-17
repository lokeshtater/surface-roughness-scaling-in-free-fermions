import numpy as np
import matplotlib.pyplot as plt
import sys

lat_sizes = [int(n) for n in sys.argv[1:3]]
system=sys.argv[int(len(lat_sizes)+1)]
a=sys.argv[int(len(lat_sizes)+2)] #hopping parameter
state=sys.argv[int(len(lat_sizes)+3)]

print(f"{lat_sizes} {system} {a} {state}" )


#importing data
w_sq = {}  #for raw SR data
times = {}

for N in lat_sizes:
    data = np.genfromtxt(f"./data_SR/SR_{system}_{state}_{N}_sites.csv", delimiter=',', dtype=complex, skip_header=1)
    times[N] = [t.real for t in data[:,0]]
    w_sq[N] = data[:,1]


print("Data imported")
####alpha and z:
######################


alpha_try = np.linspace(0.5,1.5,100) #trial values of alpha and z
z_try = np.linspace(0.5,1.5,100)
N_ref = lat_sizes[0]
deviate = [] #array to store deviations

t_i,t_f = 10**(0),10**(2)

for alpha in alpha_try:
    for z in z_try:
        chi = 0.0
        for N in lat_sizes:
            indices = [i for i,t in enumerate(times[N]) if t_i<t<t_f]
            for i in indices:
                j = i + int(np.round(len(times[N])*z*np.log10(N/N_ref)/np.log10(times[N][-1]/times[N][0])))
                chi = chi + abs(w_sq[N_ref][i] -(w_sq[N][j])*(N/N_ref)**(-alpha))/abs(w_sq[N_ref][i])
        deviate.append(chi)

ind = deviate.index(min(deviate))

alpha = alpha_try[int(ind/len(z_try))]
z = z_try[int(ind%len(z_try))]

print("alpha and z found")
print("Creating Plots:")

######growth exponents and plotting:

beta1 = {} 
beta2 = {}

t1,t2 = 10**(-3),10**(-1) #region1 (growth1)
t3,t4 = 10**(0),10**(1) #region2  (growth2)


colors_blue = ["#d4f1f4","#75e6da","#189ab4","#05445e"]
colors_purple = ["#dddddf","#d0c4df","#dcabdf","#c792df"]
colors_orange = ["#f5bb00","#ec9f05","#d76a03","#bf3100"]
c = 0

for N in lat_sizes:
    
    #fitting in region1
    indices = [i for i,t in enumerate(times[N]) if t1<t<t2]
    log_w_sq = np.log10([w_sq[N][i] for i in indices])
    time1 = [times[N][i] for i in indices]
    log_t = np.log10(time1)
    beta1[N] = np.polyfit(log_t,log_w_sq,1)[0] #fitting log(w^2) against log(t)
    
    w_sq_fit1 = [t**beta1[N] for t in time1]
    

    #fitting in region2

    indices = [i for i,t in enumerate(times[N]) if t3<t<t4]
    log_w_sq = np.log10([w_sq[N][i] for i in indices])
    time2 = [times[N][i] for i in indices]
    log_t = np.log10(time2)
    beta2[N] = np.polyfit(log_t,log_w_sq,1)[0] #fitting log(w^2) against log(t)
    
    w_sq_fit2 = [t**beta2[N] for t in time2]
    
  
    plt.plot([t/(N**z) for t in times[N]], [w/(N**alpha) for w in w_sq[N]],marker='o', color=colors_blue[c], label=f'N={N}') #raw data
    plt.plot([t/(N**z) for t in time1],[w/(N**alpha) for w in w_sq_fit1], linestyle='--', color=colors_purple[c], label=fr"$t^{{{round(beta1[N],3)}}}$") #fit of region1 (t<<1/J)
    plt.plot([t/(N**z) for t in time2],[w/(N**alpha) for w in w_sq_fit2], linestyle='--', color=colors_orange[c], label=fr"$t^{{{round(beta2[N],3)}}}$") #fit of region2 (t~1/J)

    c = c+1 

z = np.round(z,3)
alpha = np.round(alpha,3)

plt.xscale('log',base=10)
plt.yscale('log',base=10)

#creating plot title
title="Free fermions"
if(system=="tight_binding"): title += " tight binding"
elif(system=="long_hopping"): title += " long range hopping ({a})"
else: title += " Unknown model"

if(state=="alt"): title += " Alternating"
elif(state=="stagg"): title += " Staggered"
elif(state=="dom"): title += " Domain wall"
else: title += " Unknown"
title += " Initial state"
    
plt.xlabel(fr'$Jt/L^{{z}}$   (J = nearest neighbour hopping stregth)',fontsize=14)
plt.ylabel(fr"$w^{{2}}(L,t)/L^{{\alpha}}$",fontsize=14)
plt.title(title)
plt.text(1.0,0.0,fr"z = {{{z}}}  alpha= {{{alpha}}}",ha="right", va="bottom",transform=plt.gca().transAxes)


plt.legend()
plt.show()
