# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import scipy.special as ss

def calStrZo(t, b, w, zd):
    '''
    calculate impedance of striplines

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    b : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.
    zd : TYPE
        DESCRIPTION.

    Returns
    -------
    zo : TYPE
        DESCRIPTION.

    '''
    if t == 0:
        t = b/1e6
    x = t/b
    am = 2/(1+2*x/(3*(1-x)))
    dwbta = x/np.pi/(1-x)
    dwbtb = (x/(2-x))**2+(0.0796*x/((w/b)+1.1*x))**am
    dwbt = dwbta*(1-0.5*np.log(dwbtb))
    wpbt = (w/(b-t))+dwbt
    fsla = (8/np.pi/wpbt)+np.sqrt((8/np.pi/wpbt)**2+6.27)
    
    # Wu's stripline geometric factor F
    fsl = (0.5/np.pi)*np.log(1+(4/np.pi/wpbt)*fsla)
    zosl = (zd/2)*fsl
    zo = zosl
    return zo

def bf0(ip, z, key):
    x9 = z/2
    y9 = -1*x9*x9
    if key == 1:
        y9 = x9*x9
    j9 = 1
    rj9 = j9
    if9 = 50
    t9 = rj9
    for k9 in range(1,if9+1):
        rk9 = k9
        print('rk9: ',rk9)
        rip = ip
        print('rip: ',ip)
        rj9 = rj9*y9/rk9/(rip+rk9)
        print('rj9: ',rj9)
        s9 = t9 + rj9
        print('s9: ',s9)
        if s9 == t9:
            print(t9)
            break
        t9 = s9
        print('t9: ',t9)

def bf(v, z, key):
    if key == 1:
        # modified bessel of 1st kind Iv
        return ss.iv(v,z)
    if key == 2:
        # bessel 1st kind Jv
        return ss.jv(v,z)
    
# Constants
# CGS units 

# c - speed of light cm/s
# gam - gyromagnetic radio MHz/(100(Teslas)) 
# pi
# uo - permeability of free space H/cm
# eo - permitivity of free space F/cm
c = 2.9978e10
gam = 2.8e6
pi = np.pi
uo = 1.2566e-8
eo = 8.854e-14

## ENTER THESE VALUES ##
# 1 - stipline
# 2 - microstrip 
kstmi = 1

# number of greens function terms to use
ngft = 55

# fpms - 4pims, 4(pi)(Ms) in grams
# ef - ferrite dielectric constant
# tand - ferrite loss tangent
fpms = 1780
ef = 14
tand = 0.001

# ed - dielectric constant outside ferrite
ed = 10
########################

fm = gam*fpms
zd = 120*pi/np.sqrt(ed)

## ENTER THESE VALUES ##
# r - radius in cm
# w - line width in cm
# b - dielectric thickness in cm
# t - line thickness in cm
r = 0.3
w = 0.417
b = 0.125
t = 0
########################



zo = calStrZo(t,b,w,zd)

si = np.arcsin(w/2/r)

freql = 2
freqh = 20
freqs =0.5 

for f in np.arange(freql, freqh+freqs, freqs):
    f = f*1e9
    
    #200
    # key = 2: Bessel J (1st kind)
    # key = 1: Bessel I (1st kind modified)
    key = 2
    efc = ef-1j*ef*tand
    if (np.abs(f-fm) < 0.001e9):
        f = 1.001*fm
    ak = -1*fm/f
    u = 1.0
    ueff = (u*u-ak*ak)/u
    #print(ueff)
    if (f < 0.999*fm):
        ueff = -1*ueff
        key = 1
        efc = np.conjugate(efc)
    #print(ueff)
    zeff = np.sqrt(uo*ueff/eo/efc)
    
    sic1 = np.pi/np.sqrt(3)/1.84
    sic = sic1*zd*np.abs(ak/u)/zeff
    
    zibcwu = 2*zd/(1+(sic/si)*(sic/si))
    zibc = zibcwu*zo/zd
    
    sr = 2*np.pi*f*np.sqrt(ueff*efc)*r/c
    
    b0 = bf(0, sr, key)
    a0= -1*bf(1, sr, key)
    if key == 1:
        a0 = -1*a0
    
    c1a = si*b0/2/a0
    c1c = np.pi*zd/2/zeff
    
    al11 = 0+1j
    c11 = c1a + al11*c1c
    if f < fm:
        c11 = -1*c1a + al11*c1c
    
    c1b = 0
    c2br = 0
    c2bi = 0
    for n in range(1,ngft+1):
        rn = n
        bn = bf(n,sr,key)
        an = (sr*bf(n-1, sr, key) - rn*bf(n, sr, key))/sr
        ck1 = np.sin(rn*si)**2/(rn*rn*si)
        ck2 = rn*ak/u/sr
        ck22 = ck2*ck2
        c1bd = an*an- ck22*bn*bn
        c1b = c1b + ck1*an*bn/c1bd
        ck3 = np.cos(2*rn*np.pi/3)
        ck4 = np.sin(2*rn*np.pi/3)
        c2br = c2br + ck1*(an*bn*ck3)/c1bd
        c2bi = c2bi + ck1*ck2*bn*bn*ck4/c1bd
        
    if f < fm:
        c1b = -1*c1b
    
    c1 = c11+c1b
    c2 = c1a+c2br-al11*c2bi
    c3 = c1a+c2br+al11*c2bi    
    
    if f < fm:
        c2 = -1*c2
        c3 = -1*c3
    
    c4 = c1*c1-c2*c3
    c5 = c1*c1*c1 + c2*c2*c2 + c3*c3*c3 - 3*c1*c2*c3
    c6 = c2*c2 - c1*c3
    c7 = c3*c3 - c1*c2
    sm = np.pi*zd/zeff/c5
    al1 = 0-1j
    al2 = sm*c4
    al = 1 + al1*al2
    be = al1*sm*c6
    ga = al1*sm*c7
        
    zinwu = -1*zd+al1*2*zeff*c5/np.pi/c4
    zin = zinwu*zo/zd
    
    rin = np.real(zin)
    xin = np.imag(zin)
    
    # Salay/Peppiatt reactance sign correction
    xin = -1*xin
    
    dbrtnl = 20*np.log10(np.abs(al))
    dbisol = 20*np.log10(np.abs(be))
    dbinsl = 20*np.log10(np.abs(ga))
    
    yin = 1/zin
    g_in = np.real(yin)
    b_in = np.imag(yin)
    
    # Salay/Peppiatt admittance sign correction
    b_in = -1*b_in
    
    print('{:.2f}'.format(f/1e9), 
          '{:.2f}'.format(rin), 
          '{:.2f}'.format(xin), 
          '{:.2f}'.format(dbrtnl), 
          '{:.2f}'.format(dbisol), 
          '{:.2f}'.format(dbinsl),
          '{:.4f}'.format(g_in), 
          '{:.4f}'.format(b_in))