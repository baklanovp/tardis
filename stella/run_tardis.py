#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# %matplotlib inline
# from matplotlib import rcParams
# from matplotlib.lines import lineStyles
# import matplotlib.pyplot as plt

# if True:
#     plt.style.use('./mtl2smv2.stl')



# In[ ]:


# We import the necessary packages
import numpy as np
from tardis.io.config_reader import Configuration
from tardis.model import Radial1DModel
from tardis.io.atom_data.util import download_atom_data
from astropy import units as u
import matplotlib.pyplot as plt
import copy


# # Run TARDIS for STELLA data

# In[ ]:


# fconfig = 'tardis_stella.yml'
fconfig = 'tardis_stella_lbR500M20Ni01b6Z002ncVIUbox171E06.yml'


# In[ ]:


# !wget -q -nc https://raw.githubusercontent.com/baklanovp/tardis/tree/master/stella/tardis_stella.yml


# In[ ]:


from tardis import run_tardis


# In[ ]:


print("Run with "+fconfig)

sim = run_tardis(fconfig, 
                 virtual_packet_logging=False,
                 show_convergence_plots=False,
                 export_convergence_plots=False,
                 log_level="INFO",
                 show_progress_bars=False)  # WARNING  ERROR
# sim = run_tardis("tardis_stella.yml", 
#                  virtual_packet_logging=True,
#                  show_convergence_plots=True,
#                  export_convergence_plots=True,
#                  log_level="INFO")  # WARNING  ERROR


# In[ ]:


# pause


# ## Plotting the Spectrum
# 
# Finally, plot the generated spectrum with `matplotlib`.

# In[ ]:


import matplotlib.pyplot as plt


# In[ ]:


spectrum = sim.runner.spectrum
spectrum_virtual = sim.runner.spectrum_virtual
spectrum_integrated = sim.runner.spectrum_integrated


# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
plt.figure(figsize=(10, 6.5))

spectrum.plot(label="Normal packets")
spectrum_virtual.plot(label="Virtual packets")
spectrum_integrated.plot(label='Formal integral')

plt.xlim(500, 9000)
# plt.ylim(1e37, 0.8e39)
plt.ylim(0., 1.e38)
plt.title("TARDIS example model spectrum")
plt.xlabel("Wavelength [A]")
plt.ylabel("Luminosity density [erg/s/A]")
plt.legend()
plt.show()


# In[ ]:


spectrum.plot(label="Normal packets")


# In[ ]:


# dir(spectrum_integrated)


# In[ ]:





# In[ ]:





# In[ ]:


spectrum_integrated.plot()


# #### Compute magnitudes from TARDIS spectrum

# In[ ]:


import sys, os
sys.path.append(os.path.expanduser('/home/bakl/Sn/Release/python/pystella'))
import pystella as ps                
              
    
d = 10 # pc    
flux = spectrum_integrated.luminosity_density_nu #/ (4*np.pi*ps.phys.pc2cm(d)**2)    
sp_tardis = ps.Spectrum('tardis', freq=spectrum_integrated.frequency, flux=flux)
print(len(sp_tardis.Wl))
print(sp_tardis.Wl)
print(sp_tardis.Flux)
# Signature: sp.to_mag(b, z=0.0, d=3.08568047916849e+19, magnification=1.0)
mag_tardis = {}
for bn in ['U','B','V','R','I']:
    m = sp_tardis.to_mag(ps.band.band_by_name(bn))
    mag_tardis[bn] = m
    print(f"{bn} = {m}")


# In[ ]:


plt.plot(sp_tardis.Wl2angs, sp_tardis.FluxWl2angs)


# In[ ]:


#### Load UBVRI from tt-file Stella


# #### The STELLA VS TARDIS SEDs

# In[ ]:


# 'tardis_stella_lbR500M20Ni01b6Z002ncVIUbox171E06.yml'
mdl_name = fconfig.replace("tardis_stella_","").replace(".yml","") 
print(f"{fconfig=}  {mdl_name=}")


ph = ps.Stella(name=mdl_name).get_ph()

# ph_mags = ph.flux_to_curves(['U','B','V','R','I'])
# get line with time
for i, t in enumerate(ph.Time):
    tu = t * u.day
    if tu >= time_explosion:
        break
        
print(f"t= {ph.Time[i]} prev t= {ph.Time[i-1]}")
t, sp = ph.get_tspec(i)


# In[ ]:


fig, ax = plt.subplots(figsize=(18,8))

# ax = ps.curves_plot(ph_mags)
ax.plot(sp.Wl2angs, sp.FluxWl2angs, label='ph')
ax.plot(sp_tardis.Wl2angs, sp_tardis.FluxWl2angs, label=sp_tardis.Name)

ax.legend()
ax.set_xlabel('Wavelength [A]')
ax.set_xlim(1e2,1.5e4)
ax.set_ylabel('Flux')
ax.legend(loc=4);
ax.set_title(f"t_exp= {time_explosion} v_in= {config.model.structure.v_inner_boundary}");
ax.text(0.5, 0.97, mdl_name, horizontalalignment='center',verticalalignment='center', transform=ax.transAxes);


# #### The STELLA VS TARDIS photometry

# In[ ]:


tt = ps.Stella(name=mdl_name).get_tt()
ttinfo = tt.Info
ttdata = tt.load()
print(f"{tt.name}:  R= {ttinfo.R}  M= {ttinfo.M}  E= {ttinfo.E} (10^50 erg)")


# In[ ]:


config = Configuration.from_yaml(fconfig)
time_explosion = config.supernova.time_explosion
print(f"{time_explosion=}")

# get line with time
for i, l in enumerate(ttdata):
    t = l[0]* u.day
    if t >= time_explosion:
        break

p = ttdata[i-1]
print(f"prev t= {p[0]}")
l = ttdata[i]
mag_tt = {'U':l[7], 'B': l[8], 'V': l[9], 'I': l[10], 'R': l[11]}
print(t, mag_tt)


# In[ ]:


fig, ax = plt.subplots()

bnames = list(mag_tardis.keys())

# The width of the bars 
width = 0.15
x = np.arange(len(bnames))
performance = [10,8,6,4,2,1]

# plt.bar(y_pos, performance, align='center', alpha=0.5)

# Create the bar charts!
m_tds = [v for bv, v in mag_tardis.items()]
m_tt = [v for bv, v in mag_tt.items()]

ax.bar(x - width/2, m_tds, width, label='TARDIS')
ax.bar(x + width/2, m_tt, width, label='TT')

ax.set_ylim(-13,-17)

# ax.invert_yaxis()

ax.set_xticks(x)    # This ensures we have one tick per year, otherwise we get fewer
ax.set_xticklabels(bnames);
ax.set_ylabel('Magnitude')
ax.legend(loc=4);
ax.set_title(f"t_exp= {time_explosion} v_in= {config.model.structure.v_inner_boundary}");
ax.text(0.5, 0.97, mdl_name, horizontalalignment='center',verticalalignment='center', transform=ax.transAxes);


# In[ ]:


# import h5py

# def printall(name, obj):
#     print(name, dict(obj.attrs))


# path_mdl = './'    
# fname = os.path.expanduser( os.path.join(path_mdl, mdl_name+'.h5') )

# print(f'{fname=}')
    
# with h5py.File(fname,'r') as hf:
#     hf.visititems(printall)


# In[ ]:




