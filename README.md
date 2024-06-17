# Surface Roughness Scaling in Free Fermions

## How to use the repository:

1. Import all the scripts in a single directory.
2. Create an empty directory "data_SR" in the above directory.
3. Edit the following input parameters in main.sh:
    - list of lattice sizes
    - system type (tight binding or long range hopping)
    - initial state
    - time period of evaluation of dynamics
4. Also edit the "lattice sizes" input line in FV_exponents.py depending on the number of lattice sizes.
5. Compile and run main.sh

### Note the following:

1. The scripts are valid for free fermions. The hamiltonians available are the tight binding and long range hopping models. In case of long range hopping model, the hopping parameter $\alpha$ can be set in main.sh
2. Scripts for three initial pure states- Alternating, Staggered and Domain Wall state is provided. Modify correlators_direct.py to include more initial states.
3. Surface roughness is evaluated for half the system,i.e., bipartite particle number fluctuations are computed, directly from 2 point correlators, from the equation:
$w^2(L,t) = \sum_{i \in M}D_{ii}(t) - \sum_{i,j \in M}|D_{ij}(t)|^2$. (Where $D_{ij}(t)$ are the 2 point correlators, $M$ is half the subsystem of lattice of size $L$) This equation is only valid for fermions in a quadratic Hamiltonian. 
  
  
