"""
Python package for a 3D non-linear force-free magnetic field extrapolation.
The method: Optimization method of Wheatland et al. (2000)
The codes are referred to and edited from SSW package, which
can give preliminary nonlinfff

contact: lius@nao.cas.cn

input:
     bx by bz 3D potential field or linear force-free field
option:

output:
    bx by bz 3D nonlinear force-free field
"""
import pyfits as py
import numpy as np
import matplotlib.pyplot as plt


filebxyz=r'linearfff.fits'
hdu=py.open(filebxyz)
data=hdu[0].data
bx=data[0,:,:,:]
by=data[1,:,:,:]
bz=data[2,:,:,:]

# bx=bx.transpose(2,1,0)
# by=by.transpose(2,1,0)
# bz=bz.transpose(2,1,0)


###############################################################
# Initialization condition
bc_flag=0
potlin=1
iterations = 10000
dt = np.float64(0.00001)

afd = 1.0e-6
qsphere = 0
abs_frac_diff=1.0e-6

rsize,tsize,vsize=bx.shape
max_rst = np.max([rsize, tsize, vsize])
default_dx = 2.0/(max_rst-1)
#####################################################################
dtype=bx.dtype
x = np.empty(rsize, dtype=dtype)
y = np.empty(tsize, dtype=dtype)
z = np.empty(vsize, dtype=dtype)

offset = -default_dx*float(rsize/2)
x[:]= default_dx*np.arange(rsize)+offset
offset = -default_dx*float(tsize/2)
y[:]= default_dx*np.arange(tsize)+offset
z[:]= default_dx*np.arange(vsize)
rsize1 = rsize-1
tsize1 = tsize-1
vsize1 = vsize-1

dx=np.roll(x,-1)-x
dx[rsize-1]=dx[rsize-2]
dy=np.roll(y,-1)-y
dy[tsize-1]=dy[tsize-2]
dz=np.roll(z,-1)-z
dz[vsize-1]=dz[vsize-2]
wf=bx*0.+1.
#####################################################################



def deriv_dx(x,f):
    n=(x.shape)[0]
    n2=n-2
    nf=f.ndim
    nx, ny, nz = f.shape
    if nf !=3:
        print('F Array must be 3-d')
        return -1
    if (n < 3):
        print('Parameters must have at least 3 points')
        return -1
    if (n != nx):
        print('Mismatch between x and f[*,0,0]')
        return -1
    x12 = x-np.roll(x, -1)  #x0 - x1
    x01 = np.roll(x, 1)- x
    x02= np.roll(x, 1)-np.roll(x, -1)
    cx10 = (x12 / (x01 * x02))
    cx20 = (1. / x12 - 1. / x01)
    cx30 = (x01 / (x02 * x12))
    sx1 = (x01[1] + x02[1]) / (x01[1] * x02[1])
    sx2 = x02[1] / (x01[1] * x12[1])
    sx3 = x01[1] / (x02[1] * x12[1])
    rx1 = x12[n2] / (x01[n2] * x02[n2])
    rx2 = x02[n2] / (x01[n2] * x12[n2])
    rx3 = (x02[n2] + x12[n2]) / (x02[n2] * x12[n2])
    cx1 = np.empty((nx, ny, nz), dtype=dtype)
    cx2 = np.empty((nx, ny, nz), dtype=dtype)
    cx3 = np.empty((nx, ny, nz), dtype=dtype)
    for j in range(ny):
        for k in range(nz):
            cx1[:, j, k] = cx10
            cx2[:, j, k] = cx20
            cx3[:, j, k] = cx30
    d=f*0.
    d = np.roll(f, 1,axis=0) * cx1+f * cx2 - np.roll(f, -1, axis=0) * cx3
    d[0, :, :] = f[0, :, :] *sx1 - f[1, :, :] *sx2 + f[2, :, :] *sx3
    d[n - 1, :, :] = -f[n - 3, :, :] *rx1 + f[n - 2, :, :] *rx2 - f[n - 1, :, :] *rx3
    return d
def deriv_dy(y,f):
    n = (y.shape)[0]
    n2 = n - 2
    nf = f.ndim
    nx, ny, nz = f.shape
    if nf != 3:
        print('F Array must be 3-d')
        return -1
    if (n < 3):
        print('Parameters must have at least 3 points')
        return -1
    if (n != ny):
        print('Mismatch between x and f[0,:,0]')
        return -1
    y12 = y-np.roll(y, -1)  #x0 - x1
    y01 = np.roll(y, 1)- y
    y02= np.roll(y, 1)-np.roll(y, -1)
    cy10 = (y12 / (y01 * y02))
    cy20 = (1. / y12 - 1. / y01)
    cy30 = (y01 / (y02 * y12))
    sy1 = (y01[1] + y02[1]) / (y01[1] * y02[1])
    sy2 = y02[1] / (y01[1] * y12[1])
    sy3 = (y01[1] / (y02[1] * y12[1]))
    ry1 = y12[n2] / (y01[n2] * y02[n2])
    ry2 = y02[n2] / (y01[n2] * y12[n2])
    ry3 = (y02[n2] + y12[n2]) / (y02[n2] * y12[n2])
    cy1 = np.empty((ny, nx, nz), dtype=dtype)
    cy2 = np.empty((ny, nx, nz), dtype=dtype)
    cy3 = np.empty((ny, nx, nz), dtype=dtype)
    for j in range(nx):
        for k in range(nz):
            cy1[:, j, k] = cy10
            cy2[:, j, k] = cy20
            cy3[:, j, k] = cy30
    cy1 = cy1.transpose(1, 0, 2)
    cy2 = cy2.transpose(1, 0, 2)
    cy3 = cy3.transpose(1, 0, 2)
    d = f * 0.
    d = np.roll(f, 1, axis=1) * cy1 + f * cy2 - np.roll(f, -1, axis=1) * cy3
    d[:, 0, :] = f[:, 0, :] *sy1 - f[:, 1, :] *sy2 + f[:, 2, :] *sy3
    d[:, n - 1, :] = -f[:, n - 3, :] *ry1 + f[:, n - 2, :] *ry2 - f[:, n - 1, :] *ry3
    return d
def deriv_dz(z,f):
    n = (z.shape)[0]
    n2 = n - 2
    nf = f.ndim
    nx, ny, nz = f.shape
    if nf != 3:
        print('F Array must be 3-d')
        return -1
    if (n < 3):
        print('Parameters must have at least 3 points')
        return -1
    if (n != nz):
        print('Mismatch between x and f[0,0,:]')
        return -1
    z12 = z-np.roll(z, -1)  #x0 - x1
    z01 = np.roll(z, 1)- z
    z02= np.roll(z, 1)-np.roll(z, -1)
    cz10 = (z12 / (z01*z02))
    cz20 = (1./z12 - 1./z01)
    cz30 = (z01 / (z02 * z12))
    sz1 = (z01[1]+z02[1])/(z01[1]*z02[1])
    sz2 = z02[1]/(z01[1]*z12[1])
    sz3 = z01[1]/(z02[1]*z12[1])
    rz1 = z12[n2]/(z01[n2]*z02[n2])
    rz2 = z02[n2]/(z01[n2]*z12[n2])
    rz3 = (z02[n2]+z12[n2]) / (z02[n2]*z12[n2])
    cz1 = np.empty((nz, ny, nx), dtype=dtype)
    cz2 = np.empty((nz, ny, nx), dtype=dtype)
    cz3 = np.empty((nz, ny, nx), dtype=dtype)
    for j in range(ny):
        for k in range(nx):
            cz1[:, j, k] = cz10
            cz2[:, j, k] = cz20
            cz3[:, j, k] = cz30
    cz1 = cz1.transpose(2, 1, 0)
    cz2 = cz2.transpose(2, 1, 0)
    cz3 = cz3.transpose(2, 1, 0)
    d = np.roll(f, 1, axis=2) * cz1 + f * cz2 - np.roll(f, -1, axis=2) * cz3
    d[:, :, 0] = f[:, :, 0] *sz1 - f[:, :, 1] *sz2 + f[:, :, 2] *sz3
    d[:, :, n - 1] = -f[:, :, n - 3] *rz1 + f[:, :, n - 2] *rz2 - f[:, :, n - 1] *rz3
    return d

def curl_xyz(ax,ay,az,x,y,z):
    curl_x = deriv_dy(y, az) - deriv_dz(z, ay)
    curl_y = deriv_dz(z, ax) - deriv_dx(x, az)
    curl_z = deriv_dx(x, ay) - deriv_dy(y, ax)
    return curl_x,curl_y,curl_z
def div_xyz(ax,ay,az,x,y,z):
    daxdx = deriv_dx(x, ax)
    daydy = deriv_dy(y, ay)
    dazdz = deriv_dz(z, az)
    div_a = daxdx + daydy + dazdz
    return div_a
def cross_xyz(ax,ay,az,bx,by,bz):
    cx = ay * bz - az * by
    cy = az * bx - ax * bz
    cz = ax * by - ay * bx
    return cx,cy,cz

def vector_ops_fff(wf,x,y,z,bx,by,bz,):
    curl_bx,curl_by,curl_bz=curl_xyz(bx,by,bz,x,y,z)
    div_b=div_xyz(bx,by,bz,x,y,z)
    omega_x,omega_y,omega_z=cross_xyz(curl_bx,curl_by,curl_bz,bx,by,bz)
    b2 = (bx ** 2 + by ** 2 + bz ** 2)
    omega2 = (omega_x ** 2 + omega_y ** 2 + omega_z ** 2) ###add myself
    ok = np.where(b2 > 0.0)
    if ok == 0:
        print('No nonzero B')
        return curl_bx,curl_by,curl_bz,div_b,omega_x,omega_y,omega_z,b2
    omega_x[ok] = (omega_x[ok] - div_b[ok] * bx[ok]) / b2[ok]
    omega_y[ok] = (omega_y[ok] - div_b[ok] * by[ok]) / b2[ok]
    omega_z[ok] = (omega_z[ok] - div_b[ok] * bz[ok]) / b2[ok]
    omega2 = (omega_x ** 2 + omega_y ** 2 + omega_z ** 2) * wf
    omega_x = omega_x * wf
    omega_y = omega_y * wf
    omega_z = omega_z * wf
    return curl_bx,curl_by,curl_bz,div_b,omega_x,omega_y,omega_z,b2,omega2

def obj_funct_fff(bx, by, bz, rsize, tsize, vsize, rsize1, tsize1, vsize1, b2, omega2, dx, dy, dz):
    a = np.array(dx[0])
    b = dx[0:rsize - 2] + dx[1:rsize - 1]
    c = np.array(dx[rsize - 2])
    # coeffx =np.concatenate((a,b,c),axis=0)
    a = np.append(a, b)
    a = np.append(a, c)
    coeffx = a * 0.5

    a = np.array(dy[0])
    b = dy[0:tsize - 2] + dy[1:tsize - 1]
    c = np.array(dy[tsize - 2])
    # coeffx =np.concatenate((a,b,c),axis=0)
    a = np.append(a, b)
    a = np.append(a, c)
    coeffy = a * 0.5

    a = np.array(dz[0])
    b = dz[0:vsize - 2] + dz[1:vsize - 1]
    c = np.array(dz[vsize - 2])
    # coeffx =np.concatenate((a,b,c),axis=0)
    a = np.append(a, b)
    a = np.append(a, c)
    coeffz = a * 0.5

    coeffxp3 = np.empty((rsize,tsize,vsize),dtype=dtype)
    coeffyp3 = np.empty((tsize,vsize,rsize),dtype=dtype)
    coeffzp3 = np.empty((vsize,tsize,rsize),dtype=dtype)

    for i in range(vsize):
        for j in range(tsize):
            coeffxp3[:, j, i] = coeffx
    for i in range(rsize):
        for j in range(vsize):
            coeffyp3[:, j, i] = coeffy
    for i in range(rsize):
       for j in range(tsize):
           coeffzp3[:, j, i] = coeffz
    coeffxp3 = coeffxp3
    coeffyp3 = coeffyp3.transpose((2, 0, 1))
    coeffzp3 = coeffzp3.transpose((2, 1, 0))
    dV = (coeffxp3) * (coeffyp3) * (coeffzp3)
    ok = np.where(b2 > 0.0)
    if ok == 0:
        print('No nonzero B')
        return -1
    a=dV * b2 * omega2
    l= np.sum(dV * b2 * omega2)
    return l

def force_fff(omega_x, omega_y, omega_z, curl_bx, curl_by, curl_bz,div_b,omega2,b2,bx, by, bz,x,y,z):
    #First curl of omega cross b
    temp_x, temp_y, temp_z = cross_xyz(omega_x, omega_y, omega_z, bx, by, bz)
    f_x, f_y, f_z = curl_xyz(temp_x, temp_y, temp_z, x, y, z)
    #Next subtract omega cross curl B
    temp_x, temp_y, temp_z = cross_xyz(omega_x, omega_y, omega_z, curl_bx, curl_by, curl_bz)
    f_x = f_x - (temp_x)
    f_y = f_y - (temp_y)
    f_z = f_z - (temp_z)
    #Next subtract grad omega dot B
    #temp = dot_xyz(omega_x, omega_y, omega_z, bx, by, bz)
    temp = omega_x*bx+omega_y*by+omega_z*bz
    #grad_xyz, temp, x, y, z, temp_x, temp_y, temp_z
    temp_x = deriv_dx(x, temp)
    temp_y = deriv_dy(y, temp)
    temp_z = deriv_dz(z, temp)
    f_x = f_x - (temp_x)
    f_y = f_y - (temp_y)
    f_z = f_z - (temp_z)
    #Add omega times div_b
    f_x = f_x + omega_x*div_b
    f_y = f_y + omega_y*div_b
    f_z = f_z + omega_z*div_b
    #Add omega^2 times B
    f_x = f_x + omega2*bx
    f_y = f_y + omega2*by
    f_z = f_z + omega2*bz
    return f_x,f_y,f_z

def evolve_fff(i,wf,x,y,z,bx,by,bz,curl_bx, curl_by, curl_bz,div_b,omega_x, omega_y, omega_z,b2,omega2,rsize,tsize,vsize,dx, dy, dz,dt,l_value):
    afd = 1.0e-6
    rsize2 = rsize-1
    tsize2 = tsize-1
    vsize2 = vsize-1

    tbx= bx.copy()
    tby= by.copy()
    tbz= bz.copy()

    dbdt_x, dbdt_y, dbdt_z=force_fff(omega_x, omega_y, omega_z, curl_bx, curl_by, curl_bz,div_b,omega2,b2,bx, by, bz,x,y,z)
    localdt = dt #ok
    done = 0.

    while((localdt > dt/10000.0) and (done ==0)):
        bx[1:rsize2, 1:tsize2, 1:vsize2] = bx[1:rsize2, 1:tsize2, 1:vsize2] + dbdt_x[1:rsize2, 1:tsize2, 1:vsize2] * localdt
        by[1:rsize2, 1:tsize2, 1:vsize2] = by[1:rsize2, 1:tsize2, 1:vsize2] + dbdt_y[1:rsize2, 1:tsize2, 1:vsize2] * localdt
        bz[1:rsize2, 1:tsize2, 1:vsize2] = bz[1:rsize2, 1:tsize2, 1:vsize2] + dbdt_z[1:rsize2, 1:tsize2, 1:vsize2] * localdt
        curl_bx, curl_by, curl_bz, div_b, omega_x, omega_y, omega_z, b2, omega2 = vector_ops_fff(wf, x, y, z, bx, by, bz)
        local_l = obj_funct_fff(bx, by, bz, rsize, tsize, vsize, rsize1, tsize1, vsize1, b2, omega2, dx, dy, dz)

        frac_diff = np.abs((local_l - l_value) / l_value)
        if(frac_diff < afd):
            print('Converged; frac_diff = ', frac_diff)
            done = 1.
            convergence_flag = 1.
            dt = localdt
            return local_l,convergence_flag,bx,by,bz,curl_bx, curl_by, curl_bz,div_b,omega_x, omega_y, omega_z,b2,omega2,dt

        if (local_l < l_value):
            done = 1.
            localdt = localdt * 1.01
        if (local_l >= l_value):
            bx = tbx.copy()
            by = tby.copy()
            bz = tbz.copy()

            curl_bx, curl_by, curl_bz, div_b, omega_x, omega_y, omega_z, b2, omega2 = vector_ops_fff(wf, x, y, z, bx,by, bz)

            if(localdt > dt/10000.0):
                done = 0.
                localdt = localdt / 2.0
                print('Reducing dt to ' , localdt)
            if (localdt <= dt/10000.0):
                done = 1.
    if(localdt <= dt/10000.0):
        convergence_flag = 1.
        print('NOT converging')
    if (localdt > dt/10000.0):
        convergence_flag = 0.
    dt = localdt
    return local_l,convergence_flag,bx,by,bz,curl_bx, curl_by, curl_bz,div_b,omega_x, omega_y, omega_z,b2,omega2,dt

convergence_flag = 0.
for i in range(iterations):
    if (convergence_flag == 1):
        print('It seems to have converged or')
        print('is not decreasing for small time step')
        break
    if (i ==0 ):
        curl_bx,curl_by,curl_bz,div_b,omega_x,omega_y,omega_z,b2,omega2=vector_ops_fff(wf,x,y,z,bx,by,bz)
        l_value = obj_funct_fff(bx, by, bz, rsize, tsize, vsize, rsize1, tsize1, vsize1, b2, omega2, dx, dy, dz)
        delta_l = 1.0e20

    local_l, convergence_flag, bx, by, bz, curl_bx, curl_by, curl_bz, div_b, omega_x, omega_y, omega_z, b2, omega2, dt = evolve_fff(i,wf,x,y,z,bx,by,bz,curl_bx, curl_by, curl_bz,div_b,omega_x, omega_y, omega_z,b2,omega2,rsize,tsize,vsize,dx, dy, dz,dt,l_value)
    delta_l = (l_value - local_l) / l_value
    print(i, local_l, l_value)
    if (convergence_flag ==0 and local_l < l_value):
        l_value = local_l
    i_final = i
pbxyz=np.stack((bx,by,bz))
hdu = py.PrimaryHDU(pbxyz)
hdulist = py.HDUList([hdu])
hdulist.writeto(r'nonlinearfff.fits')
