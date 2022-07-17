# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 15:36:21 2019
@author: Dr. Christopher Hutchison 
Imperial College London 
2022

Version 0.1
"""
import numpy as np
def cart2pol(v):
    r = np.sqrt(np.square(v[0])+np.square(v[1])+np.square(v[2]))
    pol = np.arccos(np.absolute(v[2])/r)
    az = np.arctan2(v[1],v[0])
    return r, pol, az

def pol2cart(r, pol, az):
    x = r*np.sin(pol)*np.cos(az)
    y = r*np.sin(pol)*np.sin(az)
    z = r*np.cos(pol)
    return [x,y,z]
           
def indicatrix_radii(n,theta,phi):       
    a=n[0]
    b=n[1]
    c=n[2]        
    r = (a*b*c)/np.sqrt((b**2)*(c**2)*(np.sin(theta)**2)*(np.cos(phi)**2) + \
                        (a**2)*(c**2)*(np.sin(theta)**2)*(np.sin(phi)**2) + \
                        (a**2)*(b**2)*(np.cos(theta)**2))
    return r

def euler_rotation(v,the,phi,psi):
    rot = np.array[[np.cos(the)*np.cos(psi), np.sin(phi)*np.sin(the)*np.cos(psi)-np.cos(phi)*np.sin(psi), np.cos(phi)*np.sin(the)*np.cos(psi)+np.sin(phi)*np.sin(psi)],\
                   [np.cos(the)*np.sin(phi), np.sin(phi)*np.sin(the)*np.sin(psi)+np.cos(phi)*np.cos(psi), np.cos(phi)*np.sin(the)*np.sin(psi)-np.sin(phi)*np.cos(psi)],\
                   [np.sin(the), np.sin(phi)*np.sin(the), np.cos(phi)*np.cos(the)]]
    new_v = np.matmul(v,rot)
    return new_v    

def rotate_vec(k,vec,theta): # Rotate point <Vec> about the vector <K> by angle <theta>
    cost = np.cos(theta)    
    sint = np.sin(theta)
    u = k[0]    
    v = k[1]
    w = k[2]
    x = vec[0]
    y = vec[1]
    z = vec[2]
    
    new_vec = np.array([\
                u*(u*x+v*y+w*z)*(1-cost) + x*cost + (-1*w*y+v*z)*sint, \
                v*(u*x+v*y+w*z)*(1-cost) + y*cost +    (w*x-u*z)*sint, \
                w*(u*x+v*y+w*z)*(1-cost) + z*cost + (-1*v*x+u*y)*sint ])
#    x =    (cost+ux*ux*(1-cost))*v[0] + (ux*uy*(1-cost)-uz*sint)*v[1] + (ux*uz*(1-cost)+uy+sint)*v[2]
#    y = (ux*uy*(1-cost)+uz*sint)*v[0] +    (cost+uy*uy*(1-cost))*v[1] + (uy*uz*(1-cost)-ux+sint)*v[2]
#    z = (uz*ux*(1-cost)-uy*sint)*v[0] + (uz*uy*(1-cost)+ux*sint)*v[1] +    (cost+uz*uz*(1-cost))*v[2]
#    new_v= np.array([x,y,z])
    return new_vec

def rotate_ind(n,theta,phi,acc):
    a=n[0]
    b=n[1]
    c=n[2]

    k  = [np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]
    k /= np.linalg.norm(k)
    #create arbirarty vector perpendicular to k
    r  = np.random.randn(3)         # take a random vector
    r -= np.multiply(r.dot(k),k)    # make it orthogonal to k
    r /= np.linalg.norm(r)          # normalize it 
    s = np.cross(k,r)               # make third orthogonal vector
    s /= np.linalg.norm(s) 
    Beta1=0
    Beta2=1e10 

    for w in np.linspace(-np.pi/2,np.pi/2,int(acc)): 
#        rtemp = np.multiply(r,np.cos(w)) + np.multiply(s,np.sin(w))
#        stemp = -1*np.multiply(r,np.sin(w)) + np.multiply(s,np.cos(w))
        rtemp = rotate_vec(k,r,w)
        stemp = rotate_vec(k,s,w)
        d1    = 1/np.sqrt((rtemp[0]**2)/(a**2) + (rtemp[1]**2)/(b**2) + (rtemp[2]**2)/(c**2))
        b1temp= np.sqrt((d1*rtemp[0])**2+(d1*rtemp[1])**2+(d1*rtemp[2])**2)
        d2    = 1/np.sqrt((stemp[0]**2)/(a**2) + (stemp[1]**2)/(b**2) + (stemp[2]**2)/(c**2))
        b2temp = np.sqrt((d2*stemp[0])**2+(d2*stemp[1])**2+(d2*stemp[2])**2)
        
        #print w, rtemp, b1temp, b2temp
        if b1temp > Beta1: #Compare and and update vector, looking for maximum value of 
            Beta1 = b1temp
            r_best = rtemp
        #if b1temp < Beta2:
            Beta2 = b2temp
            s_best = stemp

    if np.dot(r_best,s_best) > 0.01 or np.dot(r_best,k) > 0.01 or np.dot(s_best,k) > 0.01: #check vectors are perpendicular
        print('Error vectors not perpendicular')
        print(np.dot(r_best,s_best), np.dot(r_best,k), np.dot(s_best,k))
    rr,rpol,raz = cart2pol(r_best)
    sr,spol,saz = cart2pol(s_best)
    
    
    if np.abs(np.cos(spol)) > np.abs(np.cos(rpol)):
        r_final = s_best/np.linalg.norm(s_best)
        s_final = r_best/np.linalg.norm(r_best)   
    else:
        r_final = r_best/np.linalg.norm(r_best)
        s_final = s_best/np.linalg.norm(s_best)
    
    
    #r_final = r_best/np.linalg.norm(r_best)
    #s_final = s_best/np.linalg.norm(s_best)
    
    #print Beta1, Beta2
    A = Beta1
    B = Beta2
    C = A-B

    return k,r_final, s_final, [A,B,C]

def Vscan(fileloc, dPhi, dTheta, dPol, dAngle):
    print('Starting E-field decomposition')
    print('Loading Xtal parameters from: '+ fileloc)
    try:
        [n, Eps, Dp_A, Abs, xtaltype] = xtalload(fileloc)
    except:
        return 0
    print(n)
    phi     = np.linspace(0,2*np.pi,dPhi)   #Azimuth
    theta   = np.linspace(0.001,  np.pi-0.001,dTheta) #Polar (don't scan 0 and pi because isotropic on axis for Uniaxial)
    pol_A   = np.linspace(0, np.pi,dPol)   #Linear polarisation angle
    
    efield  = np.empty([dTheta,dPhi,dPol,3]) # empty array for E1, E2 & E_total
    overlap = np.empty([dTheta,dPhi,dPol,2]) # empty <(mu.e)^2>
    i=0
    n_test  = np.empty([dTheta,dPhi,3])
    
    for p in range(phi.size):
        i += 1
        print('('+str(i)+'/'+str(phi.size)+')')
        for t in range(theta.size):
            first_pol     = np.array([np.sin(theta[t]+(np.pi/2))*np.cos(phi[p]),\
                                      np.sin(theta[t]+(np.pi/2))*np.sin(phi[p]),\
                                      np.cos(theta[t]+(np.pi/2))]) #creat a vector for polaristion vertically aligned 
            first_pol    /= np.linalg.norm(first_pol)
            [K,R,S,Elip]  = rotate_ind(n,theta[t],phi[p],dAngle)
            n_test[t,p,:] = [Elip[0],Elip[1],Elip[2]]
            
            for l in range(pol_A.size):
                pol             = rotate_vec(K,first_pol,pol_A[l]) #rotation polarisation
                pol            /= np.linalg.norm(pol) #Normalise 
                efield[t,p,l,0] = np.absolute(np.dot(R,pol)) #Take dot product of polarisation with R axis
                efield[t,p,l,1] = np.absolute(np.dot(S,pol)) #Take dot product of polarisation with S axis
                efield[t,p,l,2] = np.sqrt(np.square(efield[t,p,l,0])+np.square(efield[t,p,l,1])) #Total field
                overlap[t,p,l,:] = Dp_over(xtaltype,R,S,Dp_A) #find overlap dipole orientation with optical axes
                #n_test[t,p,l,0] = np.abs((R[0]*n[0]- R[1]*n[1])) + np.abs((R[0]*n[0]- R[2]*n[2])) + np.abs((R[1]*n[1]- R[2]*n[2]))
                #n_test[t,p,l,1] = np.abs((S[0]*n[0]- S[1]*n[1])) + np.abs((S[0]*n[0]- S[2]*n[2])) + np.abs((S[1]*n[1]- S[2]*n[2]))
                
    return theta, phi, pol_A, efield, overlap, n_test

def Dp_over(xtaltype,R,S,Dp_A):
    overlap = np.empty(2)
    if 'Isotropic' in str(xtaltype): #Isotropic crystal rule
        overlap = [1,1]
    elif 'Uniaxial' in str(xtaltype): # Uniaxial crystal
        rr,rpol,raz = cart2pol(R)
        sr,spol,saz = cart2pol(S)
        
        # if np.abs(rpol-np.pi/2) < 0.01: # Exception for looking down the c (O-O) axis 
        #     overlap[:] =[(0.5)*np.square(np.sin(Dp_A[0])),(0.5)*np.square(np.sin(Dp_A[0]))]
        # else:
        # overlap[:] = [np.square(np.sin(rpol))*np.square(np.cos(Dp_A[0]))+np.square(np.cos(rpol))*(0.5)*np.square(np.sin(Dp_A[0])),\
        #                   (0.5)*np.square(np.sin(Dp_A[0]))]
        overlap[:] = [np.square(np.cos(rpol))*np.square(np.cos(Dp_A[0]))+np.square(np.sin(rpol))*(0.5)*np.square(np.sin(Dp_A[0])),\
                          (0.5)*np.square(np.sin(Dp_A[0]))]

            
    elif 'Biaxial' in str(xtaltype): # Orthorhomnic crystal
        Dv = pol2cart(1,Dp_A[0],Dp_A[1]) #Normalised Dipole vector
        Dv = Dv/np.linalg.norm(Dv)
        overlap[:] = [np.square(np.dot(R,Dv)),np.square(np.dot(S,Dv))] #Overlap of Diopole vector with optical axis
        #r,pol,az = cart2pol(v)
        #overlap = np.square(np.sin(pol+Dp_A[0]))*np.square(np.cos(az+Dp_A[1]))    
        #overlap = [1,1]
    else:
        print('error crystal type not found')
        print(xtaltype)
        return 0
    return overlap    

def save_overlap(fileloc,xtalfile,overlap, efield):
    if fileloc == '':
        print('Error: No file selected.....')
        return 0
    if u'_overlap.npz' not in fileloc:
        fileloc = fileloc + '_overlap.npz'
    np.savez(fileloc,overlap=overlap,efield=efield)
    print('Saved overlap to '+fileloc)   

def calc_energy_range(energymin,energymax,dE):
    energyrange = np.logspace(np.log10(float(energymin)),np.log10(float(energymax)),float(dE))
    return energyrange
        
def calc_energy_density(energyrange,spotx,spoty):
    area = (spotx/2)*(spoty/2)*np.pi()
    energydensity = energyrange/area
    return energydensity
    
def pulsesave(fileloc,energy,spot,wavelength,pulselength):
    if fileloc == '':
        print('Error: No file selected.....')
        return 0
    if u'_pulse' not in fileloc:
        fileloc = fileloc + '_pulse'
    if u'.cfg' not in fileloc:
        fileloc = fileloc + '.cfg'
    f = open(fileloc,'w')
    f.write('Pulse parameter file \n' + \
            'version 0.1 \n' + \
            '-------------------------- \n' + \
            'Pulse energy (J) \n' + \
            energy[0] + '\t\n' + \
            'Spot size (um) \n' + \
            spot[0] + '\t' + spot[1] + '\t\n' + \
            'Wavelength (nm) \n' + \
            wavelength + '\t\n' + \
            'Pulselength (fs) \n' + \
            pulselength + '\t\n' + \
            '')
    print('Saved pulse parameters to '+fileloc)
    f.close 
    
def pulseload(fileloc):
    if fileloc == '':
        print('Error: No file selected.....')
        return None
    f = open(fileloc,'r')
    try:
        for line in f:
            if 'Pulse parameter file' in line:
                templine = next(f)
                if not 'version 0.1' in templine:
                    print('WARNING!: version mismatch of loaded file')
            if 'Pulse energy' in line:
                templine = next(f)
                energy = np.array(float(templine))
            if 'Spot size' in line:
                templine = next(f)
                spot = np.array(templine.split())
                spot = np.array([float(i) for i in spot])
            if 'Wavelength' in line:
                templine = next(f)
                wavelength = np.array(float(templine))
            if 'Pulselength' in line:
                templine = next(f)
                pulselength = np.array(float(templine))
    except:
        print('Error: Load failed')
        f.close
        return None
    f.close
    print('Loaded:'+fileloc)
    return energy, spot, wavelength, pulselength            
        
def xtalsave(fileloc,n,Eps,Dp_A,Abs,xtaltype):
    if fileloc == '':
        print('Error: No file selected.....')
        return 0
    if u'_xtal' not in fileloc:
        fileloc = fileloc + '_xtal'
    if u'.cfg' not in fileloc:
        fileloc = fileloc + '.cfg'
    f = open(fileloc,'w')
    f.write('Xtal parameter file \n' + \
            'version 0.1 \n' + \
            '-------------------------- \n' + \
            'Xtal type \n' + \
            str(xtaltype) + '\n' + \
            'Refractive index \n' + \
            str(n[0]) + '\t' + str(n[1]) + '\t' + str(n[2]) +'\t\n' + \
            'Extinction coefficients (S0, S1, S2) \n' + \
            str(Eps[0]) + '\t' + str(Eps[1]) + '\t' + str(Eps[2]) +'\t\n' + \
            'Dipole angles (Polar, Azimuth)) \n' + \
            str(Dp_A[0]) + '\t' + str(Dp_A[1]) + '\t\n' + \
            'Absorption \n' + \
            str(Abs) + '\t\n' + \
            '')
    print('Saved Xtal parameters to ' +fileloc)
    f.close   

def xtalload(fileloc):
    if fileloc == '':
        print('Error: No file selected.....')
        return None
    f = open(fileloc,'r')
    #try:
    for line in f:
        if 'Xtal parameter file' in line:
            templine = next(f)
            if not 'version 0.1' in templine:
                print('WARNING!: version mismatch of loaded file')
        if 'Xtal type' in line:
            templine = next(f)
            xtaltype = templine.split()
        if 'Refractive index' in line:
            templine = next(f)
            n = templine.split()
            n = np.array([float(i) for i in n])
        if 'Extinction coefficients' in line:
            templine = next(f)
            Eps = templine.split()
            Eps = np.array([float(i) for i in Eps])
        if 'Dipole angles' in line:
            templine = next(f)
            Dp_A = templine.split()
            Dp_A = np.array([float(i) for i in Dp_A])
        if 'Absorption' in line:
            templine = next(f)
            Abs = np.array(templine.split())
#    except:
#        print('Error: Load failed')
#        f.close
#        return None
    f.close
    print('Loaded:'+fileloc)
    return n, Eps, Dp_A, Abs, xtaltype