import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
import datetime
import sys



def collector(fname_z, fname_v, corr_thresh=60, lam_a = 1.5, lam_v = 1):
    # here get out the base timestamp from the filename
    o = fname_v[-10:-3]
    print(o)
    h = o[:2]
    m = o[2:4]
    s = o[4:]
    s = float(s) + (float(h)*60+float(m))*60

    # s is now the base time to be added to all timestamps in the data.

    # now load in the data and reshape / clean it.
    dat = np.loadtxt(fname_v)
    header = ['time','counter','status',
              'u','v','w1','w2','amp1',
              'amp2','amp3','amp4','snr1',
              'snr2','snr3','snr4','cor1',
              'cor2','cor3','cor4']
    dat = pd.DataFrame(dat,columns=header)
    dat.time += s
    dat = dat.drop(['status','counter'],axis=1)
    dat['w']=(dat.w1+dat.w2)/2.0

    # now load in the elevation timeseries and append it to the proper locations
    ele = pd.read_csv(fname_z)
    t,z = (ele.time_s, ele.position_mm)
    f = interp1d(t,z,fill_value='extrapolate')
    dat['z']=f(dat.time)
    
    # now flip dat['z']...
    dat.z = dat.z.max() - dat.z
    # and convert it to m
    dat.z = dat.z/1000

    # remove all points with insufficient correlation
    thresh = corr_thresh # correlation threshold
    n0 = dat.shape[0]    
    mask = (dat.cor1>=thresh)&(dat.cor2>=thresh)&(dat.cor3>=thresh)&(dat.cor4>=thresh)
    dat = dat[mask]
    n1 = dat.shape[0]
    print('collector filtered {}% of values for insufficient correlation'.format(round(100*(n0-n1)/n0,2)))

    # now remove points with spuriously high accelerations
    g = 9.81
    # check for spuriously high accelerations
    mask_a = ((dat[['u','v','w']].iloc[1:]/np.diff(dat.time)[:,np.newaxis]).abs() > lam_a*g).any(1)
    # now check for spuriously high velocities simultaneously
    vmag =  np.linalg.norm(dat[['u','v','w']],axis=1)[1:]
    vmean = vmag.mean()
    vstd = vmag.std()
    mask_v = vmag > vmean+lam_v*vstd
    mask1 = mask_a & mask_v 
    mask = np.ones(len(dat)).astype('bool')
    # find locations where both the velocity and acceleration are spurious
    mask[1:]=~mask1
    # filter
    n0 = dat.shape[0]    
    dat = dat[mask]
    n1 = dat.shape[0]
    print('collector filtered {}% of values for excessive acceleration'.format(round(100*(n0-n1)/n0,2)))

    # now save out the data at a new filename
    dat.to_csv(fname_v[:-4]+'-cleaned.dat',index=False)
    return dat




def binner(data,filename,nbins,scale='geometric'):

    # notice here the lower limit of hte binning is 1e-3.
    if scale=='geometric':
        zbins = np.geomspace(1e-3,data.z.max(),nbins+1) # the set of geometric bins
        digits = np.digitize(data.z,zbins[:-1]) # the digits of the z values into the geometric bins
        zbar = ([(zbins[i-1]*zbins[i])**0.5 for i in digits]) # add binned values to the dataset
        zbins = (zbins[1:-1]*zbins[:-2])**0.5 # the bins afte rcontraction with geometric means

    elif scale=='linear':
        zbins = np.linspace(1e-3,data.z.max(),nbins+1)
        digits = np.digitize(data.z,zbins[:-1])
        zbar = ([(zbins[i-1]+zbins[i])*0.5 for i in digits])
        zbins = (zbins[1:-1] + zbins[:-2])*0.5
        
    ubar = np.array([data.u[digits==i].mean() for i in range(1,nbins)]) # the binned velocities
    data['zbar']=zbar #  this is a bit janky but it seems fine.    
    cleaned_filename = filename[:-7]+'cleaned.csv'
    print('binned z values written to '+cleaned_filename) 
    # save that dataframe . . . 
    data.to_csv(cleaned_filename, index=False)
    
    # now calculate flow statistics
    # average flows vs elevation
    ubar = np.array([data.u[digits==i].mean() for i in range(1,nbins)]) # the binned velocities
    vbar = np.array([data.v[digits==i].mean() for i in range(1,nbins)]) # all components
    wbar = np.array([data.w[digits==i].mean() for i in range(1,nbins)])
    # fluctuations vs elevation
    up = np.array([data.u[digits==i].std() for i in range(1,nbins)]) # the binned rms velocities
    vp = np.array([data.v[digits==i].std() for i in range(1,nbins)]) # all components
    wp = np.array([data.w[digits==i].std() for i in range(1,nbins)])
    # the reynolds stresses u'w' (tau=-\rho u'w')
    f = lambda i: (( data.u[digits==i] - data.u[digits==i].mean() )*
                   ( data.w[digits==i]- data.w[digits==i].mean() )
                  ).mean() # utility function
    upwp = np.array([f(i) for i in range(1,nbins)])

    # turbulent kinetic energy
    tke = 0.5*(up**2+vp**2+wp**2)

    # write all values to a new dataframe
    out = np.array([zbins,ubar,vbar,wbar,up,vp,wp,tke,upwp]).T.tolist()
    header = ['z','ubar','vbar','wbar','up','vp','wp','tke','upwp']
    out = pd.DataFrame(out,columns=header)
    
    # save the new dataframe
    analysis_filename = filename[:-7]+'analyzed.csv'
    out.to_csv(analysis_filename, index=False) 
    print('flow statistics written to '+ analysis_filename)
    return out





def profile_fitter(data,stat,filename, zmin=2,zmax=5):
    """
    data -- output of the collector
    stat -- output of the binner
    filename -- used to title the plot
    zmin -- smallest z value for fitting
    zmax -- largest z value for fitting

    returns: shear velocity (cm/s) and shields number ( )
    """

    z,u = (stat.z,stat.ubar)
    z = z*100 # to cm
    u = u*100 # to cm/s
    
    # first plot the mean profile
    plt.scatter(z,u,edgecolor='darkred',facecolor='none', label='Mean Profile')
    plt.xscale('log')
    dz = (z.max()*z.min())**0.5 # units in cm - span of the plot
    plt.xlim(z.min(),z.max()) # units in cm
    plt.ylim(0,120) # units in cm/s
    
    # now plot the turbulent points
    # units in cm and cm/s
    plt.scatter(data.z*100,data.u*100,s=0.1,color='grey',zorder=-1,alpha=0.3, label='Instantaneous')
    
    # finally do the log fitting
    mask = (z>zmin)&(z<zmax)&(~np.isnan(z))&(~np.isnan(u))
    plt.axvline(zmin,color='black',lw=0.5)
    plt.axvline(zmax,color='black',lw=0.5)
    a,b = np.polyfit(np.log(z[mask]),u[mask],1)
    ustar = a*0.4 # shear velocity in cm/s. 0.4 is von karman constant
    rho,rhos,g,D = (1000,2650,9.81,5e-3) # parameters for shields. all in SI
    sh = rho*(ustar/100)**2/((rhos-rho)*g*D) # shields number. note unit conversion.
    
    # plot the log fit 
    plt.plot(z,a*np.log(z)+b,color='black',zorder=1,linestyle=':',lw=2, label='Log Fit')
    
    ustar1 = round(ustar,2)
    sh1 = round(sh,3)
    plt.text(2e-1,80,'Shear Vel: {} cm/s \nShields: {}'.format(ustar1,sh1),ha='left',fontsize=11)
    
    plt.xlabel('Height [cm]',fontsize=12)
    plt.ylabel('Velocity [cm/s]',fontsize=12)
    
    plt.title(filename)

    # control where the axis ticks show up
    xticks = [0.1, 0.5,1.0,2.0,4.0,8.0]
    xticklabels = [str(x) for x in xticks]
    ax = plt.gca()
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    
    vel_fig_name = filename[:-7]+'mean_profile.png'
    plt.savefig(vel_fig_name,dpi=200,facecolor='white')

    plt.show()


    # this plots all the flow statistics together

    fig,axes = plt.subplots(4,2,tight_layout=1,figsize=(8,8),sharex=True)
    keys = ['ubar','vbar','wbar','up','vp','wp','tke','upwp']
    labels= [r"$\bar{u}$", r"$\bar{v}$", r"$\bar{w}$", r"$u'$", r"$v'$", r"$w'$", r"$TKE = \frac{1}{2}(u'^2+v'^2+w'^2)$", r"$\tau = -\rho u' w'$"]
    z = stat.z*100 # to cm 
    for i,ax in enumerate(axes.flatten()):
        if i<7:
            ax.plot(z,stat[keys[i]],color='blue')
        else:
            rho=1000
            ax.plot(z,-rho*stat[keys[i]],color='blue')
            
        if i<3:
            ax.set_ylim(-0.1,1.0)
            ax.set_ylabel("[m/s]")
            ax.plot(z,stat[keys[i]],color='blue')
        elif (i>=3)&(i<6):
            ax.set_ylim(-0.01,0.2)
            ax.set_ylabel("[m/s]")
        elif i==6:
            ax.set_ylim(-0.01,0.2**2)
            ax.set_ylabel(r"[$m^2$/$s^2$]")
            ax.set_xlabel("Height [cm]")
            
        elif i==7:
            ax.set_ylim(-0.1,4)
            ax.set_ylabel("[Pa]")
            ax.set_xlabel("Height [cm]")
        
        ax.set_title(labels[i],fontsize=12)
        plt.suptitle("Flow Statistics from "+filename,fontsize=14)

        ax.set_xlim(1e-1,z.max())
        figname = filename[:-7]+'_flowstats.png'
        
    plt.savefig(figname,facecolor='white',dpi=200)
    plt.show()
        
    print('figure of flow statistics saved.')






    return ustar, sh


def power_spectrum(data,zbar):

    df = data.groupby(['zbar','t'],as_index=False).mean(0) # averages over repeated values of t at a given zbar

    # ok. now get the spectrum at a given value of t... 

    df = df[df.zbar==zbar].sort_values('t') # sorted on the time at this z value.

    # now get time and velocity
    t = df.t
    u = df.u

    # now get the spectrum
    from scipy.signal import lombscargle
    num_ls_freqs = 1000
    trange = t.max()-t.min()
    dt_min = t.diff().min()
    ls_min_freq = 1.0 / trange
    ls_max_freq = 1.0 / dt_min
    freqs = np.linspace(ls_min_freq, ls_max_freq, num_ls_freqs)*2*np.pi
    S = lombscargle(t, u,freqs)
    S = S/S.max()

    plt.loglog(freqs,S)
    plt.plot(freqs,freqs**(-5/3)*S.max()/freqs.min()**(-5/3),color='black')
    plt.ylabel("S [$m^2$]",fontsize=12)

    plt.xlabel("Frequency [Hz]",fontsize=12)
    
    plt.text(5e-2,1e-3,r'$z={}$ cm'.format(round(zbar,3)),fontsize=11)
    
