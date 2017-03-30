# canal

program to calculate electrical conductivity from molecular simulations using root-mean square displacments

this is similar to what has been described here:
    http://aip.scitation.org/doi/abs/10.1063/1.4890741

Those guys were, to the best of my knowledge, the first ones who did a spatial decomposition analysis of the electrical conductivity.

There is not much documentation. The directory data contains three example inputs that cover what is possible with canal at the moment

# Compilation

All necessary files are in the src subdirectory, including a makefile that might need changing depending on what compiler you prefer
For gcc you should be fine by simple typing make.

# input options

input needed:
    mandatory:
        - xcom/ycom/zcom:   x,y,z coordinate files for all fragments' centers of mass (in Angstrom)
                            the format is one line per snapshot, with values for each fragment separated by spaces
        - chgs:             charges of the fragments
                            the format is one line per fragment

    optional:
        - cell:             cell dimensions (if a spatial decomposition is what you're looking for)
        - nrestart:         number of restart points along the trajectory

        needed for correct conductivty calculation:
            - temp:         temperature (in K)
            - timestep:     timestep (in fs)
            - avvol:        average volume

        needed for splitting conductivity into parts from anion-anion, anion-cation, cation-cation:
            - split         1 number, turned on if != 0

        needed for spatial decomposition:
            - spatial:      3 numbers:
                            number of distance bins <-- if > 1 analysis turned on
                            starting distance
                            bin width
