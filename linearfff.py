"""
Python package for a 3D linear force-free magnetic field extrapolation.
The method: Green’s function calculation of Chiu & Hilton (1977)
The codes are referred to and edited from SSW package, which
can give preliminary linfff

contact: lius@nao.cas.cn


input:
    fits file of bx by bz 2D image
option:
    potlin==0 calculate the potential field with alpha=0
    potlin==1 calculate linear force free field with alpha=average(alpha on boundary)
output:
    bx by bz 3D potential field or linear force-free field
"""
import pyfits as py
import numpy as np
import matplotlib.pyplot as plt

filebx=r'bx.fits'
fileby=r'by.fits'
filebz=r'bz.fits'
#hdu=fits.open(p)
hdu=py.open(filebx)
bx0=hdu[0].data
hdu=py.open(fileby)
by0=hdu[0].data
hdu=py.open(filebz)
bz0=hdu[0].data
bx0 = bx0.transpose((1,0))### note IDL diff from Python
by0 = by0.transpose((1,0))
bz0 = bz0.transpose((1,0))


dtype=bx0.dtype
rsize,tsize=bx0.shape
vsize = np.min([rsize, tsize])
# vsize=16

bx =np.empty((rsize,tsize,vsize),dtype=dtype)*0.
by =np.empty((rsize,tsize,vsize),dtype=dtype)*0.
bz= np.empty((rsize,tsize,vsize),dtype=dtype)*0.

bx[:,:,0]=bx0
by[:,:,0]=by0
bz[:,:,0]=bz0

###############################################################
# Initialization condition
potlin=1
max_rst = np.max([rsize, tsize, vsize])
default_dx = 2.0/(max_rst-1)
#####################################################################

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

bzmax = np.max(bz[:,:,0])

if bzmax > 0:
    bx = bx / bzmax
    by = by / bzmax
    bz = bz / bzmax
alpha=0.

if potlin==1:
    bydx  = np.gradient(by[:,:,0], x, axis = 0)
    bxdy  = np.gradient(bx[:,:,0], y, axis = 1)
    tmp_bz = bz[:,:,0]
    xxx = np.where(tmp_bz == 0)
    if xxx != 0:
        tmp_bz[xxx] = 1.0e20
    loc_alph = (bydx-bxdy)/tmp_bz
    if xxx != 0:
        loc_alph[xxx] = 0.0
    sum = np.sum(np.abs(bz[:,:,0]))
    if sum == 0:
        alpha=0
    if sum != 0:
        alpha = np.sum(loc_alph*np.abs(bz[:,:,0]))/sum


print('Least-squares alpha= ',alpha)
if alpha !=0:
    print('Linear FFF implemented, alpha=',alpha)
else:
    print('Potential Field implemented')


pbx =np.empty((rsize,tsize,vsize),dtype=dtype)*0.
pby =np.empty((rsize,tsize,vsize),dtype=dtype)*0.
pbz= np.empty((rsize,tsize,vsize),dtype=dtype)*0.

coeffxp = np.empty((rsize,tsize),dtype=dtype)*0.
coeffyp = np.empty((rsize,tsize),dtype=dtype)*0.

a=np.array(dx[0])
b=dx[0:rsize - 2] + dx[1:rsize -1 ]
c=np.array(dx[rsize-2])
#coeffx =np.concatenate((a,b,c),axis=0)
a=np.append(a,b)
a=np.append(a,c)
coeffx=a*0.5

a=np.array(dy[0])
b=dy[0:tsize - 2] + dy[1:tsize -1 ]
c=np.array(dy[tsize-2])
#coeffx =np.concatenate((a,b,c),axis=0)
a=np.append(a,b)
a=np.append(a,c)
coeffy=a*0.5

coeffxp = np.empty((rsize,tsize,vsize),dtype=dtype)*0.
coeffyp = np.empty((tsize,vsize,rsize),dtype=dtype)*0.

for i in range(vsize):
    for j in range(tsize):
        coeffxp[:,j,i]=coeffx
for i in range(rsize):
    for j in range(vsize):
        coeffyp[:,j,i]=coeffy
coeffxp=coeffxp
coeffyp=coeffyp.transpose((2,0,1))

xp = np.empty((rsize,tsize,vsize),dtype=dtype)*0.
yp = np.empty((tsize,vsize,rsize),dtype=dtype)*0.
zp = np.empty((vsize,rsize,tsize),dtype=dtype)*0.
for i in range(vsize):
    for j in range(tsize):
        xp[:,j,i]=x
for i in range(rsize):
    for j in range(vsize):
        yp[:,j,i]=y
for i in range(tsize):
    for j in range(rsize):
        zp[:,j,i]=z

xp=xp
yp=yp.transpose((2,0,1))
zp=zp.transpose((1,2,0))

cos_az = np.cos(alpha * zp)
sin_az = np.sin(alpha * zp)
bzp= np.empty((rsize,tsize,vsize),dtype=dtype)
for i in range(vsize):
    bzp[:,:,i]=bz[:,:,0]

for j in range(tsize):
    for i in range(rsize):
        bigR=np.sqrt((x[i]-xp)**2+(y[j]-yp)**2)
        r = np.sqrt(bigR**2 + zp**2)
        cos_ar = np.cos(alpha * r)
        sin_ar = np.sin(alpha * r)
        xxx = np.where(bigR == 0)
        if xxx != 0:
            bigR[xxx] = 1.0e20
        yyy= np.where(r == 0)
        if yyy != 0:
            r[yyy] = 1.0e20
        gamma = zp*cos_ar/(bigR*r)-cos_az/bigR
        dgammadz = cos_ar*(1.0/(bigR*r)-zp ** 2/(bigR*r ** 3))-alpha*zp ** 2*sin_ar/(bigR*r**2)+alpha*sin_az/bigR
        gx = bzp * ((x[i] - xp) * dgammadz / bigR + alpha * gamma * (y[j] - yp) / bigR)
        gy = bzp * ((y[j] - yp) * dgammadz / bigR - alpha * gamma * (x[i] - xp) / bigR)
        gz = bzp*(zp*cos_ar/r**3+alpha*zp*sin_ar/r**2)
        a = bzp[20, 30, 3] * ((x[i] - xp[20, 30, 3]) * dgammadz[20, 30, 3] / bigR[20, 30, 3] + alpha * gamma[20, 30, 3] * (y[j] - yp[20, 30, 3]) / bigR[20, 30, 3])
        gx[i, j, :] = 0.0
        gy[i, j, :] = 0.0
        gz[i, j, :] = 0.0
        pbx[i, j, :] = np.sum(np.sum(coeffxp * coeffyp * gx,axis=1), axis=0) / np.pi/2.
        pby[i, j, :] = np.sum(np.sum(coeffxp * coeffyp * gy, axis=1), axis=0) / np.pi/2.
        pbz[i, j, :] = np.sum(np.sum(coeffxp * coeffyp * gz, axis=1), axis=0) / np.pi/2.
        if (i == 20 and j %8==0):
            print(pbx[i,j,3],pby[i,j,3],pbz[i,j,3])

bx[:, :, 1:vsize] =pbx[:, :, 1:vsize]
by[:, :, 1:vsize] =pby[:, :, 1:vsize]
bz[:, :, 1:vsize] =pbz[:, :, 1:vsize]
############################################################################
pbxyz=np.stack((bx,by,bz))
hdu = py.PrimaryHDU(pbxyz)
hdulist = py.HDUList([hdu])
hdulist.writeto(r'linearfff.fits')