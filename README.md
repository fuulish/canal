# canal

Calculates electrical conductivity from molecular simulations using root-mean square displacements and optionally splits it into contributions
from the individual species, and/or decomposes it spatially.

This is similar to what has been described here:
    http://aip.scitation.org/doi/abs/10.1063/1.4890741

Those guys were, to the best of my knowledge, the first ones who did a spatial decomposition analysis of the electrical conductivity.

There is not much documentation. The directory data contains three example inputs that cover what is possible with canal at the moment

# Compilation

All necessary files are in the src subdirectory, including a makefile that might need changing depending on what compiler you prefer
For gcc you should be fine by simple typing make.

# Input options

You need to create an input file that contains a couple of keywords and values in the style:

    keyword = value1 value2

The location of the input file has to be passed as first command-line argument.

All keywords only take one value, except for 'spatial', which takes one to three values.
The keywords are:

    mandatory:
        - xcom/ycom/zcom:   x,y,z coordinate files for all fragments' centers of mass (in Angstrom)
                            the format is one line per snapshot, with values for each fragment separated by spaces
        - chgs:             charges of the fragments
                            the format is one line per fragment

    optional:
        - cell:             cell dimensions (if a spatial decomposition is what you're looking for)
        - nrestart:         number of restart points along the trajectory
        - fitoffset:        where to start linear fit (fraction of data, default 0.1)
        - fitlength:        length of data for linear fit (fraction of data, default 0.9)
                            fitoffset + fitlength should not exceed 1.0, otherwise accessing unallocated memory

        needed for correct conductivty calculation:
            - temp:         temperature (in K)
            - timestep:     timestep (in fs)
            - avvol:        average volume (in Angstrom^3)

        needed for splitting conductivity into parts from anion-anion, anion-cation, cation-cation:
            - split         1 number, turned on if != 0

        needed for spatial decomposition:
            - spatial:      3 numbers:
                            number of distance bins <-- if > 1 analysis turned on
                            starting distance
                            bin width
