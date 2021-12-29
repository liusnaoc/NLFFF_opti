# NLFFF_opti
NLFFF_optimization\
contact: lius@nao.cas.cn
# linearff.py

Python package for a 3D linear force-free magnetic field extrapolation.\
The method: Green's function calculation of Chiu & Hilton (1977)\
The codes are referred to and edited from SSW package, which
can give preliminary linfff

input:\
    fits file of bx by bz 2D image\
option:\
    potlin==0 calculate the potential field with alpha=0\
    potlin==1 calculate linear force free field with alpha=average(alpha on boundary)\
output:\
    bx by bz 3D potential field or linear force-free field

# nonlinearff.py

Python package for a 3D non-linear force-free magnetic field extrapolation.\
The method: Optimization method of Wheatland et al. (2000)\
The codes are referred to and edited from SSW package, which
can give preliminary nonlinfff

input:\
     bx by bz 3D potential field or linear force-free field\
option:

output:\
    bx by bz 3D nonlinear force-free field

# Data for test
bx.fits\
by.fits\
bz.fits
