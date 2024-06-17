#!/bin/bash

systems=("tight_binding" "long_hopping") #available systems (all for free fermions)
states=("alt" "stagg" "dom") #available initial states- alternating,staggered and domain wall

lat_sizes=(50 75) #choose lattice sizes
t_i=-3 #dynamics would be computed from 10^{t_i} to 10^{t_f}
t_f=3

#chose system and hopping parameter
system=${systems[$1]}
state=${states[$1]}
alpha=1.5 #hopping parameter

for N in "${lat_sizes[@]}"; do
    echo "Computing dynamics for N=$N"
    python3 correlators_direct.py "$N" "$system" "$alpha" "$state" "$t_i" "$t_f"
    python3 SR_frm_corr.py "$N" "$system" "$state"
done

echo "Plotting..."
python3 FV_exponents.py "${lat_sizes[@]}" "$system" "$alpha" "$state"
