# -*- coding: utf-8 -*-
"""
@author: Dr. Christopher Hutchison 
Imperial College London 
2022

Version 0.1
"""
#import matplotlib.pyplot as plt
#import mpl_toolkits.mplot3d.axes3d as axes3d
import numpy as np
import E_functions

def overlap_test(n,d_phi,d_theta,plot=None):
    theta = np.linspace(0,np.pi,d_theta)
    phi = np.linspace(0,2*np.pi,d_phi)
    K = np.empty((d_theta,d_phi,3))
    R = np.empty((d_theta,d_phi,3))
    S = np.empty((d_theta,d_phi,3))
    Elip = np.empty((d_theta,d_phi,2))
    X = np.empty(0)
    Y = np.empty(0)
    Z = np.empty(0)
#    plt_lim = np.max(n)
#    fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
#    ax = fig.add_subplot(111, projection='3d')    
    i=0
    for t in range(theta.size):
        for p in range(phi.size):
            Rad = E_functions.indicatrix_radii(n,theta[t],phi[p])
            tempx = Rad * np.sin(theta[t]) * np.cos(phi[p])
            tempy = Rad * np.sin(theta[t]) * np.sin(phi[p])
            tempz = Rad * np.cos(theta[t])
            X = np.append(X, tempx)
            Y = np.append(Y, tempy)
            Z = np.append(Z, tempz)
            k, r, s, A, B =E_functions.rotate_ind(n,theta[t],phi[p])
            K[t,p,:] = k
            R[t,p,:] = r
            S[t,p,:] = s
            Elip[t,p,:] = [A,B]
            
            #Efield = 
#            ax.cla()
#            ax.set_xlim3d(-plt_lim, plt_lim)
#            ax.set_ylim3d(-plt_lim, plt_lim)
#            ax.set_zlim3d(-plt_lim, plt_lim)
#            ax.plot([0,tempx],[0,tempy],[0,tempz],'-b')
#            ax.plot([0,r[0]],[0,r[1]],[0,r[2]],'-r')
#            ax.plot([0,s[0]],[0,s[1]],[0,s[2]],'-g')
#            plt.show()
            i+=1
            print('('+str(i)+'/'+str(d_theta*d_phi)+')',end="\r")
#    fig2 = plt.figure(figsize=plt.figaspect(1))        
#    ax2 = fig2.add_subplot(111, projection='3d')
#    ax2.set_xlim3d(-plt_lim, plt_lim)
#    ax2.set_ylim3d(-plt_lim, plt_lim)
#    ax2.set_zlim3d(-plt_lim, plt_lim)
#    ax2.scatter(X,Y,Z)
    
    return theta, phi, K, R, S, Elip

def zscan_calc_scan(dE, dT, dz, energy, spot, wavelength, pulselength, epsilon, OD, overlap):
    
    # Constants     
    h = 6.626068e-34 #m2 kg / s
    c = 2.99792458e+8 #m/s
    NA = 6.0221415e+23 #molecules/mol
    
    PW = float(pulselength[0])*1e-15
    PW_2 = np.square(PW)
    Area=np.pi*(float(spot[0])/2*1e-4)*(float(spot[1])/2*1e-4)

    # Correct Epsilon from intergrated averaged Extinction coeff. to linear ♦
    epsilon_0 = (float(epsilon[0])/0.3)*overlap
    epsilon_1 = (float(epsilon[1])/0.3)*overlap
    epsilon_2 = (float(epsilon[2])/0.3)*overlap

    pathlength = 1

    Abs = float(OD)
    con_l = 1e-3*Abs/(epsilon_0)/pathlength# mol/cm3
    
    photon_energy = c*h/(float(wavelength[0])*1e-9)
    omega = (2*np.pi*c)/(float(wavelength[0])*1e-9)
    hb_om = (h/(2*np.pi))*omega

    print(omega)
    print(hb_om)
#    QE = 0.2; %Quantum efficiency


    energy_range = np.logspace(np.log10(float(energy[0])),np.log10(float(energy[1])),dE)
    energy_density = energy_range/Area
    esize = energy_density.size
    
    photon_density = energy_range/photon_energy/Area
    #photonpermole = photon_density/(con_l*NA*pathlength)


    sigma_0 = 1e3*np.log10(10)*epsilon_0/NA # 6cm2/molecule
    sigma_1 = 1e3*np.log10(10)*epsilon_1/NA # cm2/molecule
    sigma_2 = 1e3*np.log10(10)*epsilon_2/NA

    
    t     = np.linspace(-5*PW,5*PW,dT)
    tsize = t.size
    z     = np.linspace(0,pathlength,dz)
    zsize = z.size

    # Make empty arrays
    NS0_0 = np.zeros((esize,zsize))
    NS0   = np.zeros((esize, zsize  ,tsize))
    NS1   = np.zeros((esize, zsize  ,tsize))
    NS2   = np.zeros((esize, zsize  ,tsize))
    It    = np.zeros((esize, zsize+1,tsize))
    dIdE  = np.zeros((esize, zsize+1,tsize))
    Jo_cm2= np.zeros(esize)
    Ft    = np.zeros(tsize)
    
    print(esize,zsize,tsize)
    for i in range(esize): # Loop through energies
        Jo_cm2[i]  = np.log(2)*energy_range[i]/Area    ##ok<AGROW>
        It[i,0,:]= (Jo_cm2[i]/(np.sqrt(2*np.pi)*PW))*np.exp(-(np.square(t[:])/(PW_2*2))) #Make Gaussian (Might need rotate exp) 
        for j in range(zsize): # Loop through slice
            NS0_0[i,j] = con_l*NA*(pathlength/dz) #Make thin slice of sample
            NS0[i,j,0] = NS0_0[i,j]
            
            
            for k in range(1,tsize): # Loop through pulse
                Ft[k] = np.trapz(It[i,j,:k],x=t[:k],dx=dT) # Intergrate pulses into total flux as a function of time
                NS0[i,j,k] = NS0_0[i,j]*  (np.exp(-(sigma_0*Ft[k])/hb_om))              # NS0 population
                NS1[i,j,k] = NS0_0[i,j]*(1-np.exp(-(sigma_0*Ft[k])/hb_om))-NS2[i,j,k-1] # NS1 population
                NS2[i,j,k] = NS1[i,j,k]*(1-np.exp(-(sigma_2*Ft[k])/hb_om))              # NS2 population
            
            
            ###Linear Aborption only%%%
            #dIdE[i,j,:] = -NS0_0[i,j]*(sigma_0*It[i,j,:])
            ###inc Ground State Depleation%%%
            #dIdE[i,j,:]= -NS0_0[i,j]*(sigma_0*It[i,j,:])+sigma_0*NS1[i,j]*It[i,j,:] #GSD + Linear
            
            ###Linear + GSD + Excitated state Absorption%%%
            dIdE[i,j,:] = -NS0_0[i,j]*(sigma_0*It[i,j,:])-(sigma_1-sigma_0)*NS1[i,j]*It[i,j,:] #ESA + GSD + Linear
            
            ###Bleaching%%%
    #         NS0_0(i,j) = NS0_0(i,j)-NS2(i,j,end);%ESA + Bleach 
    # 
    #         dIdE(i,j,:) = -NS0_0(i,j).*(sigma_0.*It(i,j,:));% Linear
    #         dIdE1(i,j,:) = -NS0_0(i,j).*(sigma_0.*It1(i,j,:))-(sigma_1.*NS1(i,j,:).*It1(i,j,:))+(sigma_0.*NS0(i,j,:).*It1(i,j,:));%ESA + GSD + Linear
            
    #         figure(4)
    #         clf
    #         plot(squeeze(NS1(i,j,:)))
    #         hold on 
    #         plot(NS0_0(i,j))
    #         drawnow
            
            It[i,j+1,:] = It[i,j,:]+dIdE[i,j,:] #Reduce It by the amount absorbed dIdE1
    #         It1(i,j+1,:) = (It1(i,j,:)+dIdE1(i,j,:));
    #         
    #         for d = 2:length(t)
    #             dIdz_sum_temp(d)= trapz(t(1:d),It(i,j,1:d));
    #             dIdz_sum_temp1(d)= trapz(t(1:d),It1(i,j,1:d));
    #         end
    #         
    #         dIdz_sum(i,j) = dIdz_sum_temp(length(dIdz_sum_temp));
    #         dIdz_sum1(i,j) = dIdz_sum_temp1(length(dIdz_sum_temp1)); 
        print('('+str(i+1)+'/'+str(esize)+')',end="\r")
    
    #print(np.sum(NS0[:,:,-1],1))
    #print(con_l*NA*pathlength)
    #print(np.divide(np.sum(NS0[:,:,-1],1),(con_l*NA*pathlength)))
    
    #might need suming over all molecules ???? not sure
    return_NS0 = np.squeeze(np.divide(np.sum(NS0[:,:,-1],1),(con_l*NA*pathlength))) 
    return_NS1 = np.squeeze(np.divide(np.sum(NS1[:,:,-1],1),(con_l*NA*pathlength)))
    return_NS2 = np.squeeze(np.divide(np.sum(NS2[:,:,-1],1),(con_l*NA*pathlength)))
    return_energy_density = energy_density*10 #converts to mJ/mm2 from J/cm2
    return return_NS0, return_NS1, return_NS2, return_energy_density

def zscan_calc(energy, dT, dz, spot, wavelength, pulselength, epsilon, OD, overlap):
    #overlap = 1 ###############TESTING################################################################################
    # Constants     
    h = 6.626068e-34 #m2 kg / s
    c = 2.99792458e+8 #m/s
    NA = 6.0221415e+23 #molecules/mol
    
    PW = float(pulselength)*1e-15
    PW_2 = np.square(PW)
    Area=np.pi*(float(spot[0])/2*1e-4)*(float(spot[1])/2*1e-4)

    # Correct Epsilon from intergrated averaged Extinction coeff. to linear ♦
    epsilon_0 = (float(epsilon[0])/0.3333)*overlap
    epsilon_1 = (float(epsilon[1])/0.3333)*overlap
    epsilon_2 = (float(epsilon[2])/0.3333)*overlap

    pathlength = 1

    Abs = float(OD)
    con_l = 1e-3*Abs/(epsilon[0]/0.3333)/pathlength# mol/cm3
    
    photon_energy = c*h/(float(wavelength)*1e-9)
    omega = (2*np.pi*c)/(float(wavelength)*1e-9)
    hb_om = (h/(2*np.pi))*omega

    #print(omega)
    #print(hb_om)
#    QE = 0.2; %Quantum efficiency


    energy_density = (energy[0]+energy[1])/Area
    photon_density = energy/photon_energy/Area
    photonpermole = photon_density/(con_l*NA*pathlength)


    sigma_0 = 1e3*np.log10(10)*epsilon_0/NA # 6cm2/molecule
    sigma_1 = 1e3*np.log10(10)*epsilon_1/NA # cm2/molecule
    sigma_2 = 1e3*np.log10(10)*epsilon_2/NA

    
    t     = np.linspace(-5*PW,5*PW,dT)
    tsize = t.size
    z     = np.linspace(0,pathlength,dz)
    zsize = z.size

    # Make empty arrays
    NS0_0 = np.zeros((zsize))
    NS0   = np.zeros((zsize  ,tsize))
    NS1   = np.zeros((zsize  ,tsize))
    NS2   = np.zeros((zsize  ,tsize))
    It_e1 = np.zeros((zsize+1,tsize))
    It_e2 = np.zeros((zsize+1,tsize))
    dIdE_e1= np.zeros((zsize+1,tsize))
    dIdE_e2= np.zeros((zsize+1,tsize))
    Ft_e1 = np.zeros(tsize)
    Ft_e2 = np.zeros(tsize)
    #print(zsize,tsize)

    Jo_cm2  = np.log(2)*energy/Area    ##ok<AGROW>
    It_e1[0,:]= (Jo_cm2[0]/(np.sqrt(2*np.pi)*PW))*np.exp(-(np.square(t[:])/(PW_2*2))) #Make Gaussian (Might need rotate exp)
    It_e2[0,:]= (Jo_cm2[1]/(np.sqrt(2*np.pi)*PW))*np.exp(-(np.square(t[:])/(PW_2*2)))
    for j in range(zsize): # Loop through slice
        NS0_0[j] = con_l*NA*(pathlength/dz) #Make thin slice of sample
        NS0[j,0] = NS0_0[j]

        for k in range(1,tsize): # Loop through pulse
            Ft_e1[k] = np.trapz(It_e1[j,:k],x=t[:k],dx=dT) # Intergrate pulses into total flux as a function of time
            Ft_e2[k] = np.trapz(It_e2[j,:k],x=t[:k],dx=dT)
            NS0[j,k] = NS0_0[j]*  (np.exp(-(sigma_0[0]*Ft_e1[k]+sigma_0[1]*Ft_e2[k])/hb_om))            # NS0 population
            NS1[j,k] = NS0_0[j]*(1-np.exp(-(sigma_0[0]*Ft_e1[k]+sigma_0[1]*Ft_e2[k])/hb_om))-NS2[j,k-1] # NS1 population
            NS2[j,k] = NS1[j,k]*(1-np.exp(-(sigma_2[0]*Ft_e1[k]+sigma_2[1]*Ft_e2[k])/hb_om))            # NS2 population
        
        
        ###Linear Aborption only%%%
        #dIdE[i,j,:] = -NS0_0[i,j]*(sigma_0*It[i,j,:])
        ###inc Ground State Depleation%%%
        #dIdE[i,j,:]= -NS0_0[i,j]*(sigma_0*It[i,j,:])+sigma_0*NS1[i,j]*It[i,j,:] #GSD + Linear
        
        ###Linear + GSD + Excitated state Absorption%%%
        dIdE_e1[j,:] = -NS0_0[j]*(sigma_0[0]*It_e1[j,:])-(sigma_1[0]-sigma_0[0])*NS1[j]*It_e1[j,:] #ESA + GSD + Linear
        dIdE_e2[j,:] = -NS0_0[j]*(sigma_0[1]*It_e2[j,:])-(sigma_1[1]-sigma_0[1])*NS1[j]*It_e2[j,:] #ESA + GSD + Linear
        ###Bleaching%%%
#         NS0_0(i,j) = NS0_0(i,j)-NS2(i,j,end);%ESA + Bleach 
# 
#         dIdE(i,j,:) = -NS0_0(i,j).*(sigma_0.*It(i,j,:));% Linear
#         dIdE1(i,j,:) = -NS0_0(i,j).*(sigma_0.*It1(i,j,:))-(sigma_1.*NS1(i,j,:).*It1(i,j,:))+(sigma_0.*NS0(i,j,:).*It1(i,j,:));%ESA + GSD + Linear
        
       
        It_e1[j+1,:] = It_e1[j,:]+dIdE_e1[j,:] #Reduce It by the amount absorbed dIdE1
        It_e2[j+1,:] = It_e2[j,:]+dIdE_e2[j,:]
#         
#         for d = 2:length(t)
#             dIdz_sum_temp(d)= trapz(t(1:d),It(i,j,1:d));
#             dIdz_sum_temp1(d)= trapz(t(1:d),It1(i,j,1:d));
#         end
#         
#         dIdz_sum(i,j) = dIdz_sum_temp(length(dIdz_sum_temp));
#         dIdz_sum1(i,j) = dIdz_sum_temp1(length(dIdz_sum_temp1)); 
    #print('('+str(i+1)+'/'+str(esize)+')',end="\r")

    #print(np.sum(NS0[:,:,-1],1))
    #print(con_l*NA*pathlength)
    #print(np.divide(np.sum(NS0[:,:,-1],1),(con_l*NA*pathlength)))

    return_NS0 = np.squeeze(np.divide(np.sum(NS0[:,-1],0),(con_l*NA*pathlength))) 
    return_NS1 = np.squeeze(np.divide(np.sum(NS1[:,-1],0),(con_l*NA*pathlength)))
    return_NS2 = np.squeeze(np.divide(np.sum(NS2[:,-1],0),(con_l*NA*pathlength)))
    return_energy_density = energy_density*10 #converts to mJ/mm2 from J/cm2
    return return_NS0, return_NS1, return_NS2, return_energy_density