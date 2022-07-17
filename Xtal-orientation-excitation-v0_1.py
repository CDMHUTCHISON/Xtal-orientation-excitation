# -*- coding: utf-8 -*-
"""
@author: Dr. Christopher Hutchison 
Imperial College London 
2022

Version 0.1
"""

import tkinter as tk
import datetime, os
import zscan_func
import E_functions
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib import cm
plt.rcParams['image.cmap'] = 'jet'
import numpy as np

class Application:
    def createMainwindow(self):
        #Create mainframe
        self.MainFrame = tk.Frame(root,borderwidth=5,).grid(columnspan=2, rowspan=6)
        root.columnconfigure(1, weight=1)
        root.rowconfigure(1, weight=1)
           
        # Xtal parameter location and buttons
        self.xtalloc=tk.StringVar(self.MainFrame, value='')
        tk.Label(self.MainFrame, text='Xtal parameter file:').grid(row=0,column=0,sticky="nsew")
        self.Et_xtalloc = tk.Entry (self.MainFrame,textvariable=self.xtalloc,width=100).grid(row=0,column=1,columnspan=13,sticky="nsew")
        self.B_xtalloc = tk.Button (self.MainFrame, text="...", \
            command=lambda : self.setloc(self.xtalloc,'xtal')).grid(row=0,column=14,sticky="nsew")        
        self.B_xtal_gen = tk.Button (self.MainFrame, text="Gen.", command=self.xtal_parameters).grid(row=0,column=15,sticky="nsew")
        
        # Vector calculation buttons
        self.overlaploc=tk.StringVar(self.MainFrame, value='')
        tk.Label (self.MainFrame, text='Overlap array').grid(row=1,column=0)
        self.Et_vecloc = tk.Entry (self.MainFrame,textvariable=self.overlaploc,width=100).grid(row=1,column=1,columnspan=13,sticky="nsew")
        
        self.dTheta=tk.IntVar(self.MainFrame,value='21')
        tk.Label (self.MainFrame, text='dTheta',justify='right',width=5).grid(row=2,column=2)
        self.Et_dTheta = tk.Entry (self.MainFrame,textvariable=self.dTheta,width=5,justify='center').grid(row=2,column=3,sticky="nsew")
        self.dPhi=tk.IntVar(self.MainFrame,value='21')
        tk.Label (self.MainFrame, text='dPhi',justify='right',width=5).grid(row=2,column=4)
        self.Et_dPhi = tk.Entry (self.MainFrame,textvariable=self.dPhi,width=5,justify='center').grid(row=2,column=5,sticky="nsew")   
        self.dPol=tk.IntVar(self.MainFrame,value='21')
        tk.Label (self.MainFrame, text='dPol',justify='right',width=5).grid(row=2,column=6)
        self.Et_dAngle = tk.Entry (self.MainFrame,textvariable=self.dPol,width=5,justify='center').grid(row=2,column=7,sticky="nsew")
        self.dAngle=tk.IntVar(self.MainFrame,value='100')
        tk.Label (self.MainFrame, text='dAngle',justify='right',width=5).grid(row=2,column=8)
        self.Et_dAngle = tk.Entry (self.MainFrame,textvariable=self.dAngle,width=5,justify='center').grid(row=2,column=9,sticky="nsew")
        self.B_vecloc = tk.Button (self.MainFrame, text="...", \
            command=lambda : self.setloc_npz(self.overlaploc,'overlap')).grid(row=1,column=14,sticky="nsew")   
        self.vector_gen = tk.Button (self.MainFrame, text="Gen.", command=self.Vscan).grid(row=1,column=15,sticky="nsew")        
        
        # Pulse parameter location and buttons
        self.pulseloc=tk.StringVar(self.MainFrame, value='')
        tk.Label(self.MainFrame, text='Pulse parameter file:').grid(row=3,column=0,sticky="nsew")
        self.Et_pulseloc = tk.Entry (self.MainFrame,textvariable=self.pulseloc,width=100).grid(row=3,column=1,columnspan=13,sticky="nsew")
        self.B_pulseloc = tk.Button (self.MainFrame, text="...", \
            command=lambda : self.setloc(self.pulseloc,'pulse')).grid(row=3,column=14,sticky="nsew")
        #self.B_pulseloc_gen = tk.Button (self.MainFrame, text="Gen.", command=self.Pulse_parameters).grid(row=3,column=15,sticky="nsew")        
        
        # Scan parameters buttons
        tk.Label (self.MainFrame, text='Scan steps:').grid(row=4,column=0)
        self.dE=tk.IntVar(self.MainFrame,value='10')
        tk.Label (self.MainFrame, text='dE',justify='right',width=5).grid(row=5,column=2)
        self.Et_dE = tk.Entry (self.MainFrame,textvariable=self.dE,width=5,justify='center').grid(row=5,column=3,sticky="nsew")
        self.dt=tk.IntVar(self.MainFrame,value='10')
        tk.Label (self.MainFrame, text='dT',justify='right',width=5).grid(row=5,column=4)
        self.Et_dt = tk.Entry (self.MainFrame,textvariable=self.dt,width=5,justify='center').grid(row=5,column=5,sticky="nsew")   
        self.dz=tk.IntVar(self.MainFrame,value='5')
        tk.Label (self.MainFrame, text='dz',justify='right',width=5).grid(row=5,column=6)
        self.Et_dz = tk.Entry (self.MainFrame,textvariable=self.dz,width=5,justify='center').grid(row=5,column=7,sticky="nsew")
        
        # Start button        
        self.B_Start = tk.Button(self.MainFrame,text="Angle Scan",command = self.Zscan).grid(row=6,column=0,sticky="nsew")
        # Quit button
        self.B_Start2 = tk.Button(self.MainFrame,text="Energy Scan",command = self.Z_energy_scan).grid(row=6,column=1,sticky="nsew")
        self.B_QUIT = tk.Button(self.MainFrame,text="Quit",command = root.destroy).grid(row=6,column=14,columnspan=2,sticky="nsew")
    
    def setloc(self, var, string):
        folderloc = os.getcwd()
        self.x = tk.filedialog.askopenfilename(initialdir = folderloc,title = 'Load '+string+' config', \
            filetypes = (("config file","*_"+string+".cfg"),("all files","*.*")))
        var.set(self.x)      

    def setloc_npz(self, var, string):
        folderloc = os.getcwd()
        self.x = tk.filedialog.askopenfilename(initialdir = folderloc,title = 'Load '+string+' config', \
            filetypes = (("config file","*_"+string+".npz"),("all files","*.*")))
        var.set(self.x)  
        
    def Vscan(self):
        print(self.xtalloc.get())
        print(self.dPhi.get())
        print(self.dTheta.get())
        print(self.dPol.get())
        print(self.dAngle.get())
        try:
            theta, phi, pol_A, efield, overlap, n_test = E_functions.Vscan(self.xtalloc.get(),self.dPhi.get(), self.dTheta.get(), self.dPol.get(), self.dAngle.get())
        except:
            print('Something went wrong')
            return 0
        
        overlapsave = str(self.xtalloc.get()[:-4])+'_'+ str(self.dPhi.get()) +'_'+ str(self.dTheta.get()) + '_'+ str(self.dPol.get())+'_overlap.npz'
        self.overlaploc.set(overlapsave)
        E_functions.save_overlap(self.overlaploc.get(),self.xtalloc.get(),overlap, efield)
        
        
        [n, Eps, Dp_A, Abs, xtaltype] = E_functions.xtalload(self.xtalloc.get())
        combined  = np.multiply(np.square(efield[:,:,0,0]),overlap[:,:,0,0])
        combined1 = np.multiply(np.square(efield[:,:,0,1]),overlap[:,:,0,1])       
        
        n, Eps, Dp_A, Abs, xtaltype = E_functions.xtalload(self.xtalloc.get())
        
        fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((0, 0))
        
        
        ### Azumith vs Polarisation plot 3x3 plot ### 
        fig = plt.figure(figsize=(8,7))
        
        fig.canvas.manager.set_window_title('Azumith vs Polar (Polarisation ='+str(pol_A[0])+')')
        ax1 = fig.add_subplot(331)
        ax1.imshow(efield[:,:,0,0], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax1.title.set_text(r'$E_1$')
        ax1.set_ylabel(r'$\theta_k$')
        plt.setp(ax1,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        
        ax2 = fig.add_subplot(334)
        ax2.imshow(efield[:,:,0,1], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax2.title.set_text(r'$E_2$')
        ax2.set_ylabel(r'$\theta_k$')
        plt.setp(ax2,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        
        ax3 = fig.add_subplot(337)
        ax3.imshow(efield[:,:,0,2], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax3.title.set_text(r'$E_{total}$')
        ax3.set_ylabel(r'$\theta_k$')
        ax3.set_xlabel(r'$\phi_k$')
        plt.setp(ax3,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        
        ax4 = fig.add_subplot(332)
        ax4.imshow(overlap[:,:,0,0], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax4.title.set_text(r'$\left<\left(e_1.\mu\right)^2\right>$')
        plt.setp(ax4,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        
        ax5 = fig.add_subplot(335)
        ax5.imshow(overlap[:,:,0,1], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax5.title.set_text(r'$\left<\left(e_2.\mu\right)^2\right>$')
        plt.setp(ax5,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        ax6 = fig.add_subplot(338)
        ax6.imshow(overlap[:,:,0,0]+overlap[:,:,0,1], vmin=-0.1, vmax=1.1, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax6.title.set_text(r'$\left<\left(e_1.\mu\right)^2\right> + \left<\left(e_2.\mu\right)^2\right>$')
        ax6.set_xlabel(r'$\phi_k$') 
        plt.setp(ax6,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        ax7 = fig.add_subplot(333)
        a7 = ax7.imshow(combined, vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax7.title.set_text(r'$E_1 * \left<\left(e_1.\mu\right)^2\right>$')
        fig.colorbar(a7, ax=ax7)
        plt.setp(ax7,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        ax8 = fig.add_subplot(336)
        a8 = ax8.imshow(combined1, vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax8.title.set_text(r'$E_2 * \left<\left(e_2.\mu\right)^2\right>$')
        fig.colorbar(a8, ax=ax8)
        plt.setp(ax8,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        ax9 = fig.add_subplot(339)
        a9 = ax9.imshow(combined+combined1, vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        ax9.title.set_text(r'$E_1*\left<\left(e_1.\mu\right)^2\right> + E_2*\left<\left(e_2.\mu\right)^2\right>$')
        ax9.set_xlabel(r'$\phi_k$')
        fig.colorbar(a9, ax=ax9)
        plt.setp(ax9,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        fig.suptitle(xtaltype[0]+', '+r'$\theta_{d}$ = '+ str(Dp_A[0])+', ' +r'$\phi_{d}$ = '+ str(Dp_A[1]))
        fig.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
        
        
        # fz = 16
        # fig2 = plt.figure(figsize=(8, 6))
        # ax21 = fig2.add_subplot(1,1,1)
        # a21 = ax21.imshow(combined+combined1, vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
        # ax21.set_title(r'$E_1*\left<\left(e_1.\mu\right)^2\right> + E_2*\left<\left(e_2.\mu\right)^2\right>$',fontsize=fz)
        # ax21.set_xlabel(r'$\phi_k$',fontsize=fz)
        # ax21.set_ylabel(r'$\theta_k$',fontsize=fz)
        # fig2.colorbar(a21, ax=ax21)
        # plt.setp(ax21,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
        #          xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
        # #fig2.suptitle(xtaltype[0]+', '+r'$\theta_{d}$ = '+ str(Dp_A[0])+', ' +r'$\phi_{d}$ = '+ str(Dp_A[1]))
        # fig2.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
        
        

        
        
#        subs = np.ceil(np.sqrt(np.size(pol_A)))
#        means = np.array([])
#        fig2 = plt.figure()
#        fig3 = plt.figure()
#        for i in range(np.size(pol_A)):
#            subax = fig2.add_subplot(subs,subs,i+1)
#            subax.title.set_text(str((pol_A[i]/2/np.pi)*360))
#            combined  = np.multiply(np.square(efield[:,:,i,0]),overlap[:,:,i,0])
#            combined1 = np.multiply(np.square(efield[:,:,i,1]),overlap[:,:,i,1])
#            subax.imshow(combined, vmin=0, vmax=1, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
#            
#            subax = fig3.add_subplot(subs,subs,i+1)
#            subax.title.set_text(str((pol_A[i]/2/np.pi)*360))
#            subax.imshow(combined1, vmin=0, vmax=1, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
#            plt.show()
#            means = np.append(means,np.mean(combined+combined1))
#            print 'Average overlap :' + str(np.mean(combined+combined1))
#        print 'Total Average overlap :' + str(np.mean(means))
        
        
        ### Azumith vs Polar 3x3 plot ### 
        
        combined  = np.multiply(np.square(efield[:,0,:,0]),overlap[:,0,:,0])
        combined1 = np.multiply(np.square(efield[:,0,:,1]),overlap[:,0,:,1])
        
        fig4 = plt.figure(figsize=(8,7))
        fig4.canvas.manager.set_window_title(r'Azumith(phi) vs Polarisation(Psi) (Theta (Phi) = 0)')
        ax40 = fig4.add_subplot(3,3,1)
        ax40.imshow(efield[:,0,:,0],  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        ax40.title.set_text(r'$E_1$')
        ax40.set_ylabel(r'$\theta_k$')
        plt.setp(ax40,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        ax41 = fig4.add_subplot(3,3,2)
        ax41.imshow(overlap[:,0,:,0],  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        ax41.title.set_text(r'$\left<\left(e_1.\mu\right)^2\right>$')
        plt.setp(ax41,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        ax42 = fig4.add_subplot(3,3,3)
        a42=ax42.imshow(combined,  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        ax42.title.set_text(r'$ E_1 * \left<\left(e_1.\mu\right)^2\right>$')
        fig.colorbar(a42, ax=ax42)
        plt.setp(ax42,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        ax43 = fig4.add_subplot(3,3,4)
        ax43.imshow(efield[:,0,:,1],  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        ax43.title.set_text(r'$E_2$')
        ax43.set_ylabel(r'$\theta_k$')
        plt.setp(ax43,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        ax44 = fig4.add_subplot(3,3,5)
        ax44.imshow(overlap[:,0,:,1],  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        ax44.title.set_text(r'$\left<\left(e_2.\mu\right)^2\right>$')
        plt.setp(ax44,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        ax45 = fig4.add_subplot(3,3,6)
        a45=ax45.imshow(combined1,  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        ax45.title.set_text(r'$ E_2 * \left<\left(e_2.\mu\right)^2\right>$')
        fig.colorbar(a45, ax=ax45)
        plt.setp(ax45,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        ax46 = fig4.add_subplot(3,3,7)
        ax46.imshow(efield[:,0,:,2],  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        ax46.title.set_text(r'$E_{total}$')
        ax46.set_ylabel(r'$\theta_k$')
        ax46.set_xlabel(r'$\Psi_k$')
        plt.setp(ax46,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        ax47 = fig4.add_subplot(3,3,8)
        ax47.imshow(overlap[:,0,:,0]+overlap[:,0,:,1],  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        ax47.title.set_text(r'$\left<\left(e_1.\mu\right)^2\right> + \left<\left(e_2.\mu\right)^2\right>$')
        ax47.set_xlabel(r'$\Psi_k$')
        plt.setp(ax47,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        ax48 = fig4.add_subplot(3,3,9)
        a48 = ax48.imshow(combined+combined1,  vmin=0, vmax=1.0, extent=[theta[0],theta[-1],pol_A[0],pol_A[-1]],interpolation='nearest',aspect='auto')
        fig.colorbar(a48, ax=ax48)
        ax48.title.set_text(r'$E_1*\left<\left(e_1.\mu\right)^2\right> + E_2*\left<\left(e_2.\mu\right)^2\right>$')
        ax48.set_xlabel(r'$\Psi_k$')
        plt.setp(ax48,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
        fig4.suptitle(xtaltype[0]+', '+r'$\theta_{d}$ = '+ str(Dp_A[0])+', ' +r'$\phi_{d}$ = '+ str(Dp_A[1]))
        fig4.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
        
 
         ### 3D plots of theta vs phi
        fig3 = plt.figure(figsize=(9, 4))
        
        fig3.canvas.manager.set_window_title(r'3D plots at various polarisations')
        ax = fig3.add_subplot(131, projection='3d')
        ax1 = fig3.add_subplot(132, projection='3d')
        ax2 = fig3.add_subplot(133, projection='3d')
        
        u = np.linspace(0, 2 * np.pi, self.dPhi.get())
        v = np.linspace(0, np.pi, self.dTheta.get())
        
        # create the sphere surface
        x=n[0] * np.outer(np.cos(u), np.sin(v))
        y=n[1] * np.outer(np.sin(u), np.sin(v))
        z=n[2] * np.outer(np.ones(np.size(u)), np.cos(v))

        x=1 * np.outer(np.cos(u), np.sin(v))
        y=1 * np.outer(np.sin(u), np.sin(v))
        z=1 * np.outer(np.ones(np.size(u)), np.cos(v))
        
        combined  = np.multiply(np.square((efield[:,:,0,0])),overlap[:,:,0,0])
        combined1 = np.multiply(np.square((efield[:,:,0,1])),overlap[:,:,0,1])
 
        myheatmap = np.transpose(combined+combined1)
        ax.plot_surface(x, y, z, cstride=1, rstride=1, facecolors=cm.jet(myheatmap))
        ax.set_zlim3d(-1, 1)
        ax.set_ylim3d(-1, 1)
        ax.set_xlim3d(-1, 1)
        ax.set_box_aspect([1,1,1])
        ax.set_title(r'$\Psi_k = 0$')
        combined  = np.multiply(np.square((efield[:,:,int(self.dPol.get()/4),0])),overlap[:,:,int(self.dPol.get()/4),0])
        combined1 = np.multiply(np.square((efield[:,:,int(self.dPol.get()/4),1])),overlap[:,:,int(self.dPol.get()/4),1])
        myheatmap1 = np.transpose(combined+combined1)
        ax1.plot_surface(x, y, z, cstride=1, rstride=1, facecolors=cm.jet(myheatmap1))
        ax1.set_zlim3d(-1, 1) 
        ax1.set_ylim3d(-1, 1)
        ax1.set_xlim3d(-1, 1) 
        ax1.set_box_aspect([1,1,1])
        ax1.set_title(r'$\Psi_k = \pi/4$')
        
        combined  = np.multiply(np.square((efield[:,:,int(self.dPol.get()/2),0])),overlap[:,:,int(self.dPol.get()/2),0])
        combined1 = np.multiply(np.square((efield[:,:,int(self.dPol.get()/2),1])),overlap[:,:,int(self.dPol.get()/2),1])
        myheatmap2 = np.transpose(combined+combined1)
        ax2.plot_surface(x, y, z, cstride=1, rstride=1, facecolors=cm.jet(myheatmap2))
        ax2.set_zlim3d(-1, 1) 
        ax2.set_ylim3d(-1, 1)
        ax2.set_xlim3d(-1, 1) 
        ax2.set_box_aspect([1,1,1])
        ax2.set_title(r'$\Psi_k = \pi/2$')
        
        fig3.suptitle(xtaltype[0]+', '+r'$\theta_{d}$ = '+ str(Dp_A[0])+', ' +r'$\phi_{d}$ = '+ str(Dp_A[1]))
                
        ### Birefringence  plot ####
       
        fig5 = plt.figure(figsize=(8,4))
        fig5.canvas.manager.set_window_title(r'Relative Birefringence (minima denote optical axes)')
        gs = fig5.add_gridspec(1, 5)
        ax5 = fig5.add_subplot(gs[:2])
        a5 = ax5.pcolor(phi,theta,n_test[:,:,2])
        ax5.set_ylabel(r'$\theta_k$')
        ax5.set_xlabel(r'$\phi_k$')        
        #fig5.colorbar(a5,ax=ax5,format=fmt)
        fig5.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
        plt.setp(ax5,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])    
    
        V = n_test[:,:,2]
        V = (V-V.min())/(V.max()-V.min())
        VV = np.transpose(V)
        ax51 = fig5.add_subplot(gs[2:], projection='3d')
        ax51.plot_surface(x, y, z, cstride=1, rstride=1, facecolors=cm.jet(VV))
        ax51.set_zlim3d(-1, 1) 
        ax51.set_ylim3d(-1, 1)
        ax51.set_xlim3d(-1, 1)
        ax51.set_box_aspect([1,1,1])
        fig5.colorbar(a5,ax=ax51,format=fmt)
        fig5.suptitle(r' ')
        fig5.tight_layout(pad=0.5, w_pad=0.5, h_pad=1)
        frame1 = plt.gca()
        frame1.axes.xaxis.set_ticklabels([])
        frame1.axes.yaxis.set_ticklabels([])
        frame1.axes.zaxis.set_ticklabels([])
        
        plt.show()
        
               
    def Zscan(self):
        print('Starting scan')
        print('Loading Xtal parameters from: ' + self.xtalloc.get())
        [n, Eps, Dp_A, Abs, xtaltype] = E_functions.xtalload(self.xtalloc.get())      
        
        print('Loading Overlap array from: ' + self.overlaploc.get())
        data = np.load(self.overlaploc.get())
        overlap = data['overlap']
        efield = data['efield']
        
        print('Loading pulse parameters from: ' + self.pulseloc.get())
        [energy, spot, wavelength, pulselength] = E_functions.pulseload(self.pulseloc.get())
        #(dE, dT, dz, energy, spot, wavelength, pulselength, epsilon, OD, overlap)
        
        dimen = np.shape(overlap)
        dimax = dimen[0]*dimen[1]*dimen[2]
        NS0 = np.zeros([dimen[0],dimen[1],dimen[2],3])
        NS1 = np.zeros([dimen[0],dimen[1],dimen[2],3])
        NS2 = np.zeros([dimen[0],dimen[1],dimen[2],3])
        for i in range(dimen[0]):
            for j in range(dimen[1]):
                for k in range(dimen[2]):
                    tempE1 = np.array([np.multiply(efield[i,j,k,0],energy),0]) #Scale pulse energy by component in e1
                    output = zscan_func.zscan_calc(tempE1,self.dt.get(),self.dz.get(),spot,wavelength,pulselength,Eps,Abs,overlap[i,j,k,:])
                    NS0[i,j,k,0] = output[0]
                    NS1[i,j,k,0] = output[1]
                    NS2[i,j,k,0] = output[2]
                    
                    tempE2 = np.array([0,np.multiply(efield[i,j,k,1],energy)]) #Scale pulse energy by component in e1
                    output = zscan_func.zscan_calc(tempE2,self.dt.get(),self.dz.get(),spot,wavelength,pulselength,Eps,Abs,overlap[i,j,k,:])
                    NS0[i,j,k,1] = output[0]
                    NS1[i,j,k,1] = output[1]
                    NS2[i,j,k,1] = output[2]
                    
                    tempE3 = np.array([np.multiply(efield[i,j,k,0],energy),np.multiply(efield[i,j,k,1],energy)]) #Scale pulse energy by component in e1
                    output = zscan_func.zscan_calc(tempE3,self.dt.get(),self.dz.get(),spot,wavelength,pulselength,Eps,Abs,overlap[i,j,k,:])
                    NS0[i,j,k,2] = output[0]
                    NS1[i,j,k,2] = output[1]
                    NS2[i,j,k,2] = output[2]
            print(i*dimen[1]*dimen[2]+j*dimen[2]+k+1,'/',dimax)
        
        
        
        
        theta = np.linspace(0,np.pi,dimen[0]) 
        phi = np.linspace(0,2*np.pi,dimen[1])
        psi = np.linspace(0,np.pi,dimen[2])
        
        idx=[int(0),int(np.round(dimen[2]/4,0)),int(np.round(dimen[2]/2,0))]
        
        for indx in idx:
            fig1 = plt.figure(figsize=(10,7))
            ax11  = fig1.add_subplot(3,4,1)
            ax11.imshow(efield[:,:,indx,0], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax11.title.set_text(r'$E_1$')
            ax11.set_ylabel(r'$\theta_k$')
            plt.setp(ax11,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax12  = fig1.add_subplot(3,4,2)
            ax12.imshow(NS0[:,:,indx,0], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax12.title.set_text(r'$S_0$')
            plt.setp(ax12,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax13  = fig1.add_subplot(3,4,3)
            ax13.imshow(NS1[:,:,indx,0], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax13.title.set_text(r'$S_1$')
            plt.setp(ax13,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax14  = fig1.add_subplot(3,4,4)
            a14 = ax14.imshow(NS2[:,:,indx,0], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax14.title.set_text(r'$S_2$')
            plt.setp(ax14,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            fig1.colorbar(a14, ax=ax14)
            ax15  = fig1.add_subplot(3,4,5)
            ax15.imshow(efield[:,:,indx,1], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax15.title.set_text(r'$E_2$')
            ax15.set_ylabel(r'$\theta_k$')
            plt.setp(ax15,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax16  = fig1.add_subplot(3,4,6)
            ax16.imshow(NS0[:,:,indx,1], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax16.title.set_text(r'$S_0$')
            plt.setp(ax16,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax17  = fig1.add_subplot(3,4,7)
            ax17.imshow(NS1[:,:,indx,1], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax17.title.set_text(r'$S_1$')
            plt.setp(ax17,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax18  = fig1.add_subplot(3,4,8)
            a18 = ax18.imshow(NS2[:,:,indx,1], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax18.title.set_text(r'$S_2$')
            plt.setp(ax18,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            fig1.colorbar(a18, ax=ax18)
            ax19  = fig1.add_subplot(3,4,9)
            ax19.imshow(efield[:,:,indx,2], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax19.title.set_text(r'$E_{total}$')
            ax19.set_ylabel(r'$\theta_k$')
            ax19.set_xlabel(r'$\phi_k$')
            plt.setp(ax19,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax110  = fig1.add_subplot(3,4,10)
            ax110.imshow(NS0[:,:,indx,2], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax110.title.set_text(r'$S_0$')
            ax110.set_xlabel(r'$\phi_k$')
            plt.setp(ax110,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax111  = fig1.add_subplot(3,4,11)
            ax111.imshow(NS1[:,:,indx,2], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax111.title.set_text(r'$S_1$')
            ax111.set_xlabel(r'$\phi_k$')
            plt.setp(ax111,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            ax112  = fig1.add_subplot(3,4,12)
            a112 = ax112.imshow(NS2[:,:,indx,2], vmin=0, vmax=1.0, extent=[phi[0],phi[-1],theta[0],theta[-1]],interpolation='nearest',aspect='auto')
            ax112.title.set_text(r'$S_2$')
            ax112.set_xlabel(r'$\phi_k$')
            plt.setp(ax112,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                     xticks=[0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],xticklabels=['0','',r'$\pi$','',r'2$\pi$'])
            fig1.colorbar(a112, ax=ax112)
            fig1.suptitle(r'Populations modelling: $\theta_k$ vs $\phi_k$ ($\psi_k$ =' + '{:.2f}'.format(psi[indx])+ ')\n' + xtaltype[0]+', '+r'$\theta_{d}$ = '+ str(Dp_A[0])+', ' +r'$\phi_{d}$ = '+ str(Dp_A[1]))
            fig1.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
             
        idx=[int(0),int(np.round(dimen[1]/4,0)),int(np.round(dimen[1]/2,0))]
        
        for indx in idx:
        
            fig2 = plt.figure(figsize=(10,7))
            ax21  = fig2.add_subplot(3,4,1)
            ax21.pcolor(psi,theta,efield[:,indx,:,0],vmin=0, vmax=1.0)
            ax21.set_ylabel(r'$\theta_k$')
            ax21.title.set_text(r'$E_1$')
            plt.setp(ax21,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax22  = fig2.add_subplot(3,4,2)
            ax22.pcolor(psi,theta,NS0[:,indx,:,0],vmin=0, vmax=1.0)
            ax22.title.set_text(r'$S_0$')
            plt.setp(ax22,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax23  = fig2.add_subplot(3,4,3)
            ax23.pcolor(psi,theta,NS1[:,indx,:,0],vmin=0, vmax=1.0)
            ax23.title.set_text(r'$S_1$')
            plt.setp(ax23,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax24  = fig2.add_subplot(3,4,4)
            a24 = ax24.pcolor(psi,theta,NS2[:,indx,:,0],vmin=0, vmax=1.0)
            ax24.title.set_text(r'$S_2$')
            fig2.colorbar(a24, ax=ax24)
            plt.setp(ax24,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax25  = fig2.add_subplot(3,4,5)
            ax25.pcolor(psi,theta,efield[:,indx,:,1],vmin=0, vmax=1.0)                
            ax25.set_ylabel(r'$\theta_k$')
            ax25.title.set_text(r'$E_2$')
            plt.setp(ax25,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax26  = fig2.add_subplot(3,4,6)
            ax26.pcolor(psi,theta,NS0[:,indx,:,1],vmin=0, vmax=1.0)
            ax26.title.set_text(r'$S_0$')
            plt.setp(ax25,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax27  = fig2.add_subplot(3,4,7)
            ax27.pcolor(psi,theta,NS1[:,indx,:,1],vmin=0, vmax=1.0)
            ax27.title.set_text(r'$S_1$')
            plt.setp(ax25,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax28  = fig2.add_subplot(3,4,8)
            a28 = ax28.pcolor(psi,theta,NS2[:,indx,:,1],vmin=0, vmax=1.0)
            ax28.title.set_text(r'$S_2$')
            fig2.colorbar(a28, ax=ax28)
            plt.setp(ax25,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax29  = fig2.add_subplot(3,4,9)
            ax29.pcolor(psi,theta,efield[:,indx,:,0]+efield[:,indx,:,1],vmin=0, vmax=1.0)
            ax29.title.set_text(r'$E_{total}$')
            ax29.set_ylabel(r'$\theta_k$')
            ax29.set_xlabel(r'$\psi_k$')
            plt.setp(ax25,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax210  = fig2.add_subplot(3,4,10)
            ax210.pcolor(psi,theta,(NS0[:,indx,:,2]),vmin=0, vmax=1.0)
            ax210.title.set_text(r'$S_0$')
            ax210.set_xlabel(r'$\psi_k$')
            plt.setp(ax25,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax211  = fig2.add_subplot(3,4,11)
            ax211.pcolor(psi,theta,(NS1[:,indx,:,2]),vmin=0, vmax=1.0)
            ax211.title.set_text(r'$S_1$')
            ax211.set_xlabel(r'$\psi_k$')
            plt.setp(ax25,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            ax212  = fig2.add_subplot(3,4,12)
            a212 = ax212.pcolor(psi,theta,(NS2[:,indx,:,2]),vmin=0, vmax=1.0)
            ax212.title.set_text(r'$S_2$')
            ax212.set_xlabel(r'$\psi_k$')
            plt.setp(ax25,yticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],yticklabels=['0','',r'$\pi$/2','',r'$\pi$'],\
                 xticks=[0, np.pi/4, 2*np.pi/4, 3*np.pi/4, np.pi],xticklabels=['0','',r'$\pi$/2','',r'$\pi$'])
            fig2.colorbar(a212, ax=ax212)
            fig2.suptitle(r'Populations modelling: $\theta_k$ vs $\psi_k$ ($\phi_k$ =' + '{:.2f}'.format(phi[indx])+ ')\n' + xtaltype[0]+', '+r'$\theta_{d}$ = '+ str(Dp_A[0])+', ' +r'$\phi_{d}$ = '+ str(Dp_A[1]))
            fig2.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5) 
            
            
        plt.show()
                        
            
            
    def Z_energy_scan(self):
        print('Starting scan')
        print('Loading Xtal parameters from: ' + self.xtalloc.get())
        [n, Eps, Dp_A, Abs, xtaltype] = E_functions.xtalload(self.xtalloc.get())      
        
        print('Loading Overlap array from: ' + self.overlaploc.get())
        data = np.load(self.overlaploc.get())
        overlap = data['overlap']
        efield = data['efield']
        
        print('Loading pulse parameters from: ' + self.pulseloc.get())
        [energy, spot, wavelength, pulselength] = E_functions.pulseload(self.pulseloc.get())
        
        E_space = np.logspace(-9,-5,self.dE.get())
        NS0 = np.zeros([np.shape(E_space)[0],3])
        NS1 = np.zeros([np.shape(E_space)[0],3])
        NS2 = np.zeros([np.shape(E_space)[0],3])
        E_density = np.zeros(np.shape(E_space))
        t = 1
        p = 3
        ps = 5
        
        for i in range(np.shape(E_space)[0]):
            tempE1 = np.array([np.multiply(efield[t,p,ps,0],E_space[i]),0]) #E1 components 
            output = zscan_func.zscan_calc(tempE1, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,0] = output[0]
            NS1[i,0] = output[1]
            NS2[i,0] = output[2]
            tempE2 = np.array([0,np.multiply(efield[t,p,ps,1],E_space[i])]) #E2 components 
            output = zscan_func.zscan_calc(tempE2, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,1] = output[0]
            NS1[i,1] = output[1]
            NS2[i,1] = output[2]
            tempE3 = np.array([np.multiply(efield[t,p,ps,0],E_space[i]),np.multiply(efield[t,p,ps,1],E_space[i])])  #E_total
            output = zscan_func.zscan_calc(tempE3, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,2] = output[0]
            NS1[i,2] = output[1]
            NS2[i,2] = output[2]
            E_density[i] = output[3]
            print('('+str(i+1)+'/'+str(np.shape(E_space)[0])+')')
        
        
        fig1 = plt.figure(figsize=(12,6))
        
        ax1  = fig1.add_subplot(1,3,1)
        ax1.semilogx(E_density,NS0[:,0],'k--',E_density,NS1[:,0],'r--',E_density,NS2[:,0],'b--',\
                     E_density,NS0[:,1],'k-.',E_density,NS1[:,1],'r-.',E_density,NS2[:,1],'b-.',\
                     E_density,NS0[:,2],'k-',E_density,NS1[:,2],'r-',E_density,NS2[:,2],'b-')
        plt.legend([r'$E_1 N_{SO}$','$E_1 N_{S1}$','$E_1 N_{S2}$',\
                     '$E_2 N_{SO}$','$E_2 N_{S1}$','$E_2 N_{S2}$',\
                     '$E_{total} N_{SO}$','$E_{total} N_{S1}$','$E_{total} N_{S2}$'])
        combined  = np.multiply(np.square(efield[t,p,ps,0]),overlap[t,p,ps,0]) + np.multiply(np.square(efield[t,p,ps,1]),overlap[t,p,ps,1])
        ax1.title.set_text(r'$\left<\left(e.\mu\right)^2\right> =$ %.2f' % combined)
        ax1.set_xlabel(r'Energy density (mJ/mm$^2$)')    
        ax1.set_ylabel(r'Population')
        
        t = 3
        p = 3
        ps = 0
        
        for i in range(np.shape(E_space)[0]):
            tempE1 = np.array([np.multiply(efield[t,p,ps,0],E_space[i]),0]) #E1 components 
            output = zscan_func.zscan_calc(tempE1, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,0] = output[0]
            NS1[i,0] = output[1]
            NS2[i,0] = output[2]
            tempE2 = np.array([0,np.multiply(efield[t,p,ps,1],E_space[i])]) #E2 components 
            output = zscan_func.zscan_calc(tempE2, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,1] = output[0]
            NS1[i,1] = output[1]
            NS2[i,1] = output[2]
            tempE3 = np.array([np.multiply(efield[t,p,ps,0],E_space[i]),np.multiply(efield[t,p,ps,1],E_space[i])])  #E_total
            output = zscan_func.zscan_calc(tempE3, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,2] = output[0]
            NS1[i,2] = output[1]
            NS2[i,2] = output[2]
            E_density[i] = output[3]
            print('('+str(i+1)+'/'+str(np.shape(E_space)[0])+')')

        ax2  = fig1.add_subplot(1,3,2)
        ax2.semilogx(E_density,NS0[:,0],'k--',E_density,NS1[:,0],'r--',E_density,NS2[:,0],'b--',\
                     E_density,NS0[:,1],'k-.',E_density,NS1[:,1],'r-.',E_density,NS2[:,1],'b-.',\
                     E_density,NS0[:,2],'k-',E_density,NS1[:,2],'r-',E_density,NS2[:,2],'b-')
        # plt.legend([r'$E_1 N_{SO}$','$E_1 N_{S1}$','$E_1 N_{S2}$',\
        #              '$E_2 N_{SO}$','$E_2 N_{S1}$','$E_1 N_{S2}$',\
        #              '$E_{total} N_{SO}$','$E_{total} N_{S1}$','$E_{total} N_{S2}$'])
        combined  = np.multiply(np.square(efield[t,p,ps,0]),overlap[t,p,ps,0]) + np.multiply(np.square(efield[t,p,ps,1]),overlap[t,p,ps,1])
        ax2.title.set_text(r'$\left<\left(e.\mu\right)^2\right> =$ %.2f' % combined)
        ax2.set_xlabel(r'Energy density (mJ/mm$^2$)')
        ax2.set_ylabel(r'Population')
        
        t = 15
        p = 3
        ps = 0
        
        for i in range(np.shape(E_space)[0]):
            tempE1 = np.array([np.multiply(efield[t,p,ps,0],E_space[i]),0]) #E1 components 
            output = zscan_func.zscan_calc(tempE1, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,0] = output[0]
            NS1[i,0] = output[1]
            NS2[i,0] = output[2]
            tempE2 = np.array([0,np.multiply(efield[t,p,ps,1],E_space[i])]) #E2 components 
            output = zscan_func.zscan_calc(tempE2, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,1] = output[0]
            NS1[i,1] = output[1]
            NS2[i,1] = output[2]
            tempE3 = np.array([np.multiply(efield[t,p,ps,0],E_space[i]),np.multiply(efield[t,p,ps,1],E_space[i])])  #E_total
            output = zscan_func.zscan_calc(tempE3, self.dt.get(), self.dz.get(), spot, wavelength, pulselength, Eps, Abs, overlap[t,p,ps,:])
            NS0[i,2] = output[0]
            NS1[i,2] = output[1]
            NS2[i,2] = output[2]
            E_density[i] = output[3]
            print('('+str(i+1)+'/'+str(np.shape(E_space)[0])+')')

        ax3  = fig1.add_subplot(1,3,3)
        ax3.semilogx(E_density,NS0[:,0],'k--',E_density,NS1[:,0],'r--',E_density,NS2[:,0],'b--',\
                     E_density,NS0[:,1],'k-.',E_density,NS1[:,1],'r-.',E_density,NS2[:,1],'b-.',\
                     E_density,NS0[:,2],'k-',E_density,NS1[:,2],'r-',E_density,NS2[:,2],'b-')
        # plt.legend([r'$E_1 N_{SO}$','$E_1 N_{S1}$','$E_1 N_{S2}$',\
        #              '$E_2 N_{SO}$','$E_2 N_{S1}$','$E_1 N_{S2}$',\
        #              '$E_{total} N_{SO}$','$E_{total} N_{S1}$','$E_{total} N_{S2}$'])
        combined  = np.multiply(np.square(efield[t,p,ps,0]),overlap[t,p,ps,0]) + np.multiply(np.square(efield[t,p,ps,1]),overlap[t,p,ps,1])
        ax3.title.set_text(r'$\left<\left(e.\mu\right)^2\right> =$ %.2f' % combined)
        ax3.set_xlabel(r'Energy density (mJ/mm$^2$)')
        ax3.set_ylabel(r'Population')
        fig1.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
        plt.show()
        
    def Pulse_parameters(self):
        self.PulseWindow = tk.Toplevel(master=None)
        self.PulseWindow.title('Zscan parameter')
        
        #Energy range
        tk.Label (self.PulseWindow, text='Energy').grid(row=0,column=0)
        self.energymin=tk.StringVar(self.PulseWindow, value='1e-7')
        self.energymax=tk.StringVar(self.PulseWindow, value='1e-4')
        self.Et_energymin = tk.Entry (self.PulseWindow,textvariable=self.energymin,width=7,justify='center').grid(row=0,column=1,sticky="nsew")
        tk.Label (self.PulseWindow, text='-').grid(row=0,column=2)
        self.Et_energymax = tk.Entry (self.PulseWindow,textvariable=self.energymax,width=7,justify='center').grid(row=0,column=3,sticky="nsew")
        tk.Label (self.PulseWindow, text='J').grid(row=0,column=4)        
        
        
        #Spot size        
        self.spotx=tk.StringVar(self.PulseWindow, value='100')
        self.spoty=tk.StringVar(self.PulseWindow, value='100')
        tk.Label (self.PulseWindow, text='Spot size (FWHM)').grid(row=2,column=0)
        self.Et_spotx = tk.Entry (self.PulseWindow,textvariable=self.spotx,width=7,justify='center').grid(row=2,column=1,sticky="nsew")
        tk.Label (self.PulseWindow, text='x').grid(row=2,column=2)
        self.Et_spoty = tk.Entry (self.PulseWindow,textvariable=self.spoty,width=7,justify='center').grid(row=2,column=3,sticky="nsew")        
        tk.Label (self.PulseWindow, text='um').grid(row=2,column=4)
        
        #Area
        self.area=tk.StringVar(self.PulseWindow)
        tk.Label (self.PulseWindow, text='Area',justify='center').grid(row=3,column=0)
        self.Et_area = tk.Label (self.PulseWindow,textvariable=self.area,width=7,justify='center').grid(row=3,column=3,sticky="nsew")        
        tk.Label (self.PulseWindow, text='mm2').grid(row=3,column=4)
 
        def calc_area():
            self.area.set("%.2g" % ((float(self.spotx.get())/2000)*(float(self.spoty.get())/2000)*np.pi))
        
        #Energy density
        self.EnergyDmin=tk.StringVar(self.PulseWindow)
        self.EnergyDmax=tk.StringVar(self.PulseWindow)
        tk.Label (self.PulseWindow, text='Energy density').grid(row=4,column=0)
        self.Et_EnergyDmin = tk.Label (self.PulseWindow,textvariable=self.EnergyDmin,width=7,justify='center').grid(row=4,column=1,sticky="nsew")
        tk.Label (self.PulseWindow, text='-').grid(row=4,column=2)
        self.Et_EnergyDmax = tk.Label (self.PulseWindow,textvariable=self.EnergyDmax,width=7,justify='center').grid(row=4,column=3,sticky="nsew")        
        tk.Label (self.PulseWindow, text='mJ/mm2').grid(row=4,column=4)     
        def energy_den():
            calc_area()
            self.energy_range = E_functions.calc_energy_range(self.energymin.get(),self.energymax.get(),self.dE.get())
            self.energy_density = (self.energy_range*1000)/float(self.area.get())
            self.EnergyDmin.set("%.3f" % (self.energy_density[0]))
            self.EnergyDmax.set("%.3f" % (self.energy_density[-1]))
        energy_den()
        
        #Wavelength
        self.wavelength=tk.StringVar(self.PulseWindow,value='450')
        tk.Label (self.PulseWindow, text='Wavelength',justify='center').grid(row=5,column=0)
        self.Et_wavelength = tk.Entry (self.PulseWindow,textvariable=self.wavelength,width=7,justify='center').grid(row=5,column=3,sticky="nsew")        
        tk.Label (self.PulseWindow, text='nm').grid(row=5,column=4)
        
        #Pulse length
        self.pulselength=tk.StringVar(self.PulseWindow,value='130')
        tk.Label (self.PulseWindow, text='Pulse length',justify='center').grid(row=6,column=0)
        self.Et_pulselength = tk.Entry (self.PulseWindow,textvariable=self.pulselength,width=7,justify='center').grid(row=6,column=3,sticky="nsew")        
        tk.Label (self.PulseWindow, text='fs').grid(row=6,column=4)
 
        #Update and recalc values when pulse parameters are changed.
        def update_p(*args):
            calc_area()
            energy_den()
        
        self.energymax.trace('w',update_p)
        self.energymin.trace('w',update_p)
        self.spotx.trace('w',update_p)
        self.spoty.trace('w',update_p)
        

        
        #Load parameters
        def gui_pulseload():
            folderloc = os.getcwd()
            fileloc=tk.filedialog.askopenfilename(initialdir = folderloc,title = 'Load pulse config', \
            filetypes = (("pulse config file","*_pulse.cfg"),("all files","*.*")))       
            [energy, spot, wavelength, pulselength] = E_functions.pulseload(fileloc)
            self.energymin.set(energy[0])
            self.energymax.set(energy[1])
            self.spotx.set(spot[0])
            self.spoty.set(spot[1])
            self.wavelength.set(wavelength[0])
            self.pulselength.set(pulselength[0])
            self.pulseloc.set(fileloc)
            print('Load complete')
            self.pulseLoad=tk.Button(self.PulseWindow,text='Load',command=gui_pulseload).grid(row=7,column=1,sticky="nsew")
        
        # Save Pulse parameters
        
        def gui_pulsesave():
            folderloc = os.getcwd()
            fileloc=tk.filedialog.asksaveasfilename(initialdir = folderloc,  \
                        initialfile = str(datetime.date.today())+'_pulse.cfg', \
                        title = 'Save pulse config', filetypes = \
                        (("pulse config file","*pulse.cfg"),("all files","*.*")))
                
            E_functions.pulsesave(fileloc, \
                                [self.energymin.get(),self.energymax.get()], \
                                [self.spotx.get(),self.spoty.get()], \
                                self.wavelength.get(), \
                                self.pulselength.get())  
            self.PulseWindow.destroy()
            self.pulseloc.set(fileloc)
        self.pulseSave=tk.Button(self.PulseWindow,text='Save',command=gui_pulsesave).grid(row=7,column=3,sticky="nsew")
    
    def xtal_parameters(self):
        self.xtalWindow = tk.Toplevel(master=None)
        self.xtalWindow.title('Xtal parameter')
        self.xtal_type = tk.StringVar(self.xtalWindow,value='')
        self.Et_xtal_type = tk.Label(self.xtalWindow,textvariable=self.xtal_type)
        self.Op_xtal = tk.OptionMenu(self.xtalWindow,self.xtal_type,\
                            'Isotropic', \
                            'Uniaxial' , \
                            'Biaxial (Orthorhombic or Monoclinic)',\
                            ).grid(row=0,column=0,columnspan=7,sticky='nsew')
                    
        #self.SpaceGroup = tk.StringVar(self.xtalWindow,'')
        #tk.Label (self.xtalWindow, text='Spacegroup').grid(row=1,column=0,sticky='e')
        #self.Et_SpaceGroup = tk.Entry (self.xtalWindow, textvariable=self.SpaceGroup,width=10,justify='center').grid(row=1,column=1,sticky='w')
        
        # Refractive index
        self.n1 = tk.DoubleVar(self.xtalWindow,'1.0')
        tk.Label (self.xtalWindow, text=' n1 ').grid(row=3,column=0,sticky='e')
        self.Et_n1 = tk.Entry (self.xtalWindow, textvariable=self.n1,width=10,justify='center').grid(row=3,column=1,sticky='w')
        self.n2 = tk.DoubleVar(self.xtalWindow,'1.0')
        tk.Label (self.xtalWindow, text=' n2 ').grid(row=4,column=0,sticky='e')
        self.Et_n2 = tk.Entry (self.xtalWindow, textvariable=self.n2,width=10,justify='center').grid(row=4,column=1,sticky='w')
        self.n3 = tk.DoubleVar(self.xtalWindow,'1.0')
        tk.Label (self.xtalWindow, text=' n3 ').grid(row=5,column=0,sticky='e')
        self.Et_n3 = tk.Entry (self.xtalWindow, textvariable=self.n3,width=10,justify='center').grid(row=5,column=1,sticky='w')

        # Extinction coefficients
        self.Ep0 = tk.DoubleVar(self.xtalWindow,'44000')
        tk.Label (self.xtalWindow, text='Epsilon 0').grid(row=7,column=0,sticky='e')
        self.Et_Ep0 = tk.Entry (self.xtalWindow, textvariable=self.Ep0,width=10,justify='center').grid(row=7,column=1,sticky='w')
        self.Ep1 = tk.DoubleVar(self.xtalWindow,'28000')
        tk.Label (self.xtalWindow, text='Epsilon 1').grid(row=8,column=0,sticky='e')
        self.Et_Ep0 = tk.Entry (self.xtalWindow, textvariable=self.Ep1,width=10,justify='center').grid(row=8,column=1,sticky='w')
        self.Ep2 = tk.DoubleVar(self.xtalWindow,'200')
        tk.Label (self.xtalWindow, text='Epsilon 2').grid(row=9,column=0,sticky='e')
        self.Et_Ep0 = tk.Entry (self.xtalWindow, textvariable=self.Ep2,width=10,justify='center').grid(row=9,column=1,sticky='w')            
        
        # Dipole orientation
        self.Dp_pol = tk.DoubleVar(self.xtalWindow,'0')
        tk.Label (self.xtalWindow, text='Dipole Pol').grid(row=11,column=0,sticky='e')
        self.Et_Dp_pol = tk.Entry (self.xtalWindow, textvariable=self.Dp_pol,width=10,justify='center').grid(row=11,column=1,sticky='w')
        self.Dp_az = tk.DoubleVar(self.xtalWindow,'0')
        tk.Label (self.xtalWindow, text='Dipole Az').grid(row=12,column=0,sticky='e')
        self.Et_Dp_az = tk.Entry (self.xtalWindow, textvariable=self.Dp_az,width=10,justify='center').grid(row=12,column=1,sticky='w')                
        
        # Absorption
        self.OD = tk.DoubleVar(self.xtalWindow,'0.3')
        tk.Label (self.xtalWindow, text='Absorption').grid(row=15,column=0,sticky='e')       
        self.Et_OD = tk.Entry (self.xtalWindow, textvariable=self.OD,width=10,justify='center').grid(row=15,column=1,sticky='w')
        tk.Label (self.xtalWindow, text='OD').grid(row=15,column=2) 
        
        # Quantum yield
        self.QuantumY = tk.DoubleVar(self.xtalWindow,'0.2')
        tk.Label (self.xtalWindow, text='Q.Yield').grid(row=16,column=0,sticky='e')       
        self.Et_QuantumY = tk.Entry (self.xtalWindow, textvariable=self.QuantumY,width=10,justify='center').grid(row=16,column=1,sticky='w')
     
        def gui_xtalload():
            folderloc = os.getcwd()
            fileloc=tk.filedialog.askopenfilename(initialdir = folderloc,title = 'Load xtal config', \
            filetypes = (("config file","*_xtal.cfg"),("all files","*.*")))       
            [n, Eps, Dp_A, Abs, xtaltype] = E_functions.xtalload(fileloc)
            print(xtaltype)
            if xtaltype == 0:
                return None
            self.xtal_type(xtaltype[0])
            self.n1.set(n[0])
            self.n2.set(n[1])
            self.n3.set(n[2])
            self.Ep0.set(Eps[0])
            self.Ep1.set(Eps[1])
            self.Ep2.set(Eps[2])
            self.OD.set(Abs[0])
            print('Load complete')
        self.xtalLoad=tk.Button(self.xtalWindow,text='Load',command=gui_xtalload).grid(row=20,column=0,sticky="nsew")
        
        # Save Pulse parameters
        
        def gui_xtalsave():
            folderloc = os.getcwd()
            
            if self.n1.get() == self.n2.get() or self.n2.get() == self.n3.get() or self.n1.get() == self.n3.get():
                print(self.n1.get())
                print(self.n2.get())
                print(self.n3.get())
                self.xtal_type.set('Uniaxial')
            if (self.n1.get() == self.n2.get() and self.n2.get() == self.n3.get() and self.n1.get() == self.n3.get()):
                print(self.n1.get())
                print(self.n2.get())
                print(self.n3.get())
                self.xtal_type.set('Isotropic')         
                
            fileloc=tk.filedialog.asksaveasfilename(initialdir = folderloc,  \
                        initialfile = str(datetime.date.today())+'_'+self.xtal_type.get()+'_xtal.cfg', \
                        title = 'Save Xtal config', filetypes = \
                        (("config file","*_xtal.cfg"),("all files","*.*")))

            E_functions.xtalsave(fileloc, \
                                [self.n1.get(),self.n2.get(),self.n3.get()], \
                                [self.Ep0.get(),self.Ep1.get(),self.Ep2.get()], \
                                [self.Dp_pol.get(), self.Dp_az.get()], \
                                self.OD.get(), \
                                self.xtal_type.get(),
                                )  
            self.xtalWindow.destroy()
            self.xtalloc.set(fileloc)
        self.xtalSave=tk.Button(self.xtalWindow,text='Save',command=gui_xtalsave).grid(row=20,column=1,sticky="nsew")
        
    def __init__(self, master):
        self.MainFrame = tk.Frame(root,borderwidth=5,).grid(columnspan=2, rowspan=6)
        self.createMainwindow()        

if __name__ == '__main__':
    root = tk.Tk()
    root.title('Zscan calculator')
    app = Application(root)
    root.mainloop()