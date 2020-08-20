#!/usr/bin/env python
# coding: utf-8

# In[84]:


# Code made by Jalyn Krause from Lauren Laufmans's Mathmatica file 'AstroModels.nb'
# Edited by Jalyn Krause on 30 July 2019 - May 2020
# New code added on by Jalyn Krause May 2020 - 10 Aug. 2020

#######################################################################################
# THIS CODE WAS UPLOADED AFTER GENERATING PLOTS FOR ORIGINAL REGIONS 4-7 OF ZW049.057
#######################################################################################

# USES ARBITRARY LOCATION MODEL (APPROACH 1)

# Our objective is to measure the transmission of galaxy light through the dust clouds. 
# We work to do this at each position (x,y) in the dusty region of the image.
# We need to know the observed intensity, I, which we measure. 
# We also need an estimate of the unobscured intensity, I0, to derive the transmission. Tr(lambda) = I(lambda)/I0(lambda)y 
# We assume the galaxy is symmetric and that the outer dust plume only obscures half of the galaxy. 
# This amounts to saying the galaxy disk is tilted and we only observe dust in absorption from the side that is towards us.
# A similar feature exists on the back side of the galaxy but as it is behind most of the starlight, it produces little absorption.
# Lauren cut the galaxy image in half, aligned the two halves, and took the ratio to give us Tr(lambda, x, y).
# k(lambda) is known for our extinction curve

# In[85]:


#get_ipython().run_line_magic('matplotlib', 'notebook') # line commented out when ran in Brackets (external editor)
import sys
import sympy as sym
import scipy
from astropy.io import fits
from sympy import Symbol, exp, Eq, ImageSet, S
from sympy import solveset
from sympy import nonlinsolve
from scipy import stats
import numpy as np
#import matplotlib as plt
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
plt.rcParams.update({'figure.max_open_warning': 0})


# In[86]:


# Defines symbols used in script
tau_v = sym.Symbol('tau_v', real=True)
tau_i = sym.Symbol('tau_i', real=True)
tau_j = sym.Symbol('tau_j', real=True)
tau_h = sym.Symbol('tau_h', real=True)
tau_u = sym.Symbol('tau_u', real=True)
f_v = sym.Symbol('f_v', real=True)
f_i = sym.Symbol('f_i', real=True)
f_j = sym.Symbol('f_j', real=True)
wid = sym.Symbol('wid', real=True)


# In[87]:


# Sets up equns for plotting

# Values for total extinctions are taken from Savage & Mathis (1979) paper cited in Lauren's Senior Thesis
av = 3.10 # V band
ai = 1.87 # I band
aj = 0.87 # J band
ah = 0.614 # H band
au = 5.14 # U band 

# In terms of tau_v, prints tau values at bottom    
tau_i = solveset(Eq((ai / av) * tau_v, tau_i), tau_i) # tau_i = (ai / av) * tau_v
tau_j = solveset(Eq((aj / av) * tau_v, tau_j), tau_j) # tau_j = (aj / av) * tau_v
tau_h = solveset(Eq((ah / av) * tau_v, tau_h), tau_h) # tau_h = (ah / av) * tau_v
tau_u = solveset(Eq((au / av) * tau_v, tau_u), tau_u) # tau_u = (au / av) * tau_v

#print(tau_i, tau_j, tau_h, tau_u)


# In[88]:


# Sections 1A-J of Lauren's Original Code
# Diagram of sections shown in Fig. 3 of Lauren's Senior Thesis

# Values taken from ds9 data for 'Nick's Pillar' regions defined by Lauren L.
v_1 = [0.3221, 0.3345, 0.2844, 0.3298, 0.3554, 0.2876, 0.3708, 0.2942, 0.3089, 0.3433]
i_1 = [0.5603, 0.6047, 0.4897, 0.5459, 0.6287, 0.5270, 0.5826, 0.5520, 0.5572, 0.6138]
j_1 = [0.7325, 0.7639, 0.7007, 0.7358, 0.7444, 0.6948, 0.7413, 0.7496, 0.7607, 0.7985]

# Developed from equn (5) of Lauren's senior thesis
# V = f + (1-f)*e^(-tau_V)
# ==> tau_v = log( (V -f) / 1-f)
# f = fractional distance into the cloud where the dust is located
# Left side of equation is a ratio of the total light the observer sees in each band, in a fractional form compared to the total light emitted
# *Recall all tau values are in terms of tau_v*

def tau(f,V):
    xxx = (V-f) / (1-f)
    yyy = -1*np.log( xxx )
    yyy[xxx <= 0] = 5
    return yyy

f_input = np.linspace(0,1,101)

filenames = ['1A', '1B',  '1C', '1D', '1E', '1F', '1G', '1H', '1I', '1J']
    
for i in range(0, len(filenames)):

    #for i in range(0, len(v_1)):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(f_input, tau(f_input, v_1[i])/3.1*3.1, label="V="+str(v_1[i]))
    ax.plot(f_input, tau(f_input, i_1[i])/1.87*3.1, label="I="+str(i_1[i]))
    ax.plot(f_input, tau(f_input, j_1[i])/0.87*3.1, label="J="+str(j_1[i]))
    ax.set_ylim((0.,3))
    ax.set_xlim((0,1))
    ax.set_title(r'$\tau_v$ vs. Fractional Location in Galaxy (f)')
    ax.set_xlabel('f')
    ax.set_ylabel(r'$\tau_v$')
    ax.legend(loc = 'best')
    #print('V = ', v_1[i], 'I = ', i_1[i], 'J = ', j_1[i])
        
    #plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_1/'+ filenames[i] + '.png') #SAVES EACH FIG TO FOLDER LABELED BY REGION NAME
   


    # ADD ON TO FIND INTERSECTION POINT
    # (use link to reference later) https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gaussian_kde.html
    #X, Y = np.mgrid[0:1:100j, 0:3:100j]
    #positions = np.vstack([f_input.ravel(), (tau(f_input, 0.5603)/1.87*3.1).ravel()])
    #print(positions)
    #values = np.vstack([np.ones_like(f_input).ravel()])
    #kernel = stats.gaussian_kde(values)
    #Z = np.reshape(kernel(positions).T, X.shape)  
    
#fig.show() # (USE TO DISPLY PLOTS)
    
 


# In[89]:


# Sections 2A-I of Lauren's Original Code
# Diagram of sections shown in Fig. 3 of Lauren's Senior Thesis
# Values taken from DS9 data

v_2 = [0.3795, 0.3865, 0.3919, 0.3651, 0.3203, 0.3374, 0.2898, 0.3258, 0.3467]
i_2 = [0.7132, 0.7223, 0.6805, 0.6953, 0.6464, 0.6425, 0.6391, 0.6732, 0.6761]
j_2 = [0.8054, 0.8373, 0.8188, 0.8195, 0.7896, 0.8144, 0.8191, 0.8092, 0.7924]

# Developed from equn (5) of Lauren's senior thesis
# V = f + (1-f)*e^(-tau_V)
# ==> tau_v = log( (V -f) / 1-f)
# f = fractional distance into the cloud where the dust is located
# Left side of equation is a ratio of the total light the observer sees in each band, in a fractional form compared to the total light emitted
# *Recall all tau values are in terms of tau_v*

def tau(f,V):
    xxx = (V-f) / (1-f)
    yyy = -1*np.log( xxx )
    yyy[xxx <= 0] = 5
    return yyy

f_input = np.linspace(0,1,101)

filenames = ['2A', '2B',  '2C', '2D', '2E', '2F', '2G', '2H', '2I']
    
for j in range(0, len(filenames)):

    #for j in range(0, len(v_2)):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(f_input, tau(f_input, v_2[j])/3.1*3.1, label="V="+str(v_2[j]))
    ax.plot(f_input, tau(f_input, i_2[j])/1.87*3.1, label="I="+str(i_2[j]))
    ax.plot(f_input, tau(f_input, j_2[j])/0.87*3.1, label="J="+str(j_2[j]))
    ax.set_ylim((0.,3))
    ax.set_xlim((0,1))
    ax.set_title(r'$\tau_v$ vs. Fractional Location in Galaxy (f)')
    ax.set_xlabel('f')
    ax.set_ylabel(r'$\tau_v$')
    ax.legend(loc = 'best')
    #print('V = ', v_2[j], 'I = ', i_2[j], 'J = ', j_2[j])
    
    #plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_2/'+ filenames[j] + '.png') #SAVES EACH FIG TO FOLDER LABELED BY REGION NAME

#fig.show()  # (USE TO DISPLY PLOTS)

# In[90]:


# Sections 3A-J of Lauren's Original Code
# Diagram of sections shown in Fig. 3 of Lauren's Senior Thesis
# Values taken from DS9 data

v_3 = [0.2894, 0.3244, 0.3487, 0.4144, 0.4668, 0.3226, 0.2603, 0.2895, 0.2550, 0.3434]
i_3 = [0.4487, 0.5853, 0.6992, 0.8281, 0.8048, 0.4899, 0.4919, 0.5872, 0.4678, 0.5532]
j_3 = [0.6487, 0.7449, 0.8759, 0.9313, 0.8894, 0.6729, 0.6906, 0.7305, 0.7501, 0.7528]

# Developed from equn (5) of Lauren's senior thesis
# V = f + (1-f)*e^(-tau_V)
# ==> tau_v = log( (V -f) / 1-f)
# f = fractional distance into the cloud where the dust is located
# Left side of equation is a ratio of the total light the observer sees in each band, in a fractional form compared to the total light emitted
# *Recall all tau values are in terms of tau_v*

def tau(f,V):
    xxx = (V-f) / (1-f)
    yyy = -1*np.log( xxx )
    yyy[xxx <= 0] = 5
    return yyy

f_input = np.linspace(0,1,101)

filenames = ['3A', '3B',  '3C', '3D', '3E', '3F', '3G', '3H', '3I', '3J']
    
for k in range(0, len(filenames)):

    #for k in range(0, len(v_3)):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(f_input, tau(f_input, v_3[k])/3.1*3.1, label="V="+str(v_3[k]))
    ax.plot(f_input, tau(f_input, i_3[k])/1.87*3.1, label="I="+str(i_3[k]))
    ax.plot(f_input, tau(f_input, j_3[k])/0.87*3.1, label="J="+str(j_3[k]))
    ax.set_ylim((0.,3))
    ax.set_xlim((0,1))
    ax.set_title(r'$\tau_v$ vs. Fractional Location in Galaxy (f)')
    ax.set_xlabel('f')
    ax.set_ylabel(r'$\tau_v$')
    ax.legend(loc = 'best')
    #print('V = ', v_3[k], 'I = ', i_3[k], 'J = ', j_3[k])
    
    #plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_3/'+ filenames[k] + '.png') #SAVES EACH FIG TO FOLDER LABELED BY REGION NAME


#fig.show()  # (USE TO DISPLY PLOTS)


# In[91]:


# AREA OF SECTIONS SHOWN IN FIG. 3 OF LAUREN'S SENIOR THESIS

# Arcseconds per pixel
arcpix = 0.03962000086903572

# Area in square pixels of the regions; {1C, 1D, 3A, 3F, 3G, 3J} 
areasqpx = [28, 32, 28, 19, 17, 18]

lcm = []

# Width in arcseconds of regions
for i in range(0, len(areasqpx)):
    widtharc = arcpix * np.sqrt(areasqpx[i])

    # Arcseconds to Radians conversion: 1 rad = 180/(pi*3600 arcsec) = 206,265 arcseconds
    conv = (180*3600) / np.pi
    widthrad = widtharc / conv

    # From NED, D = 56 Mpc
    # Theta * D = ~L
    dist = 56*10**6  # [Mpc]
    
    lpc = widthrad * dist
    print('lpc = ', lpc)
    
    lcm = np.append(lcm, (3.085677*10**18) * lpc)

print('lcm = ', lcm)


# In[92]:


# FINDING EXTINCTION CURVE INFO
# Values used come from Draine 2001 paper cited in Lauren's senior thesis
# Values from method taken directly from Lauren's code, not needed to be recalculated
'''
Band	wavelength	Sigma(extinction)
U	0.365	7.62*10^-22
B	0.440	6.20*10^-22
V	0.550	4.66*10^-22
R	0.700	3.31*10^-22
I	0.900	2.27*10^-22
J	1.22	1.41*10^-22
H	1.63	8.99*10^-23
K	2.20	5.45*10^-23
L	3.45	2.37*10^-23
 	3.60	2.18*10^-23
'''    

data_1 = [[0.365, 7.62], [0.440, 6.20], [0.550, 4.66], [.700, 3.31], [0.9, 2.27], [1.22, 1.41], [1.63, 0.899]]

data_2 = [[0.7, 3.31], [0.9, 2.27], [1.22, 1.41], [1.63, 0.899], [2.20, 0.545], [3.45, 0.237], [3.60, 0.218]]

data_3 = [[0.274, 6.2], [0.344, 4.9], [0.4, 4.4], [0.44, 4.10], [0.55, 3.10], [0.70, 2.32]]

sigma_v = 4.60582*10**-22
sigma_i = 2.64231*10**-22
sigma_j = 1.35557*10**-22


# In[93]:


# CALCULATING n

# From equn. (9), (10), and (11) of Lauren's senior thesis
# L is the path length (effectively the diameter of the sphere of the dust feature)

# Taus & fs for regions:{1C, 1D, 3A, 3F, 3G, 3J} (choose them from the graphs)
# Taus = {1.351, 1.138, 1.935, 1.818, 1.312, 1.04};
# Similarly, f = {0.0456, 0.02616, 0.1785, 0.2044, 0.009954, 0.01319}

tau_reg = [1.351, 1.138, 1.935, 1.818, 1.312, 1.04]
f_reg = [0.0456,0.02616,0.1785,0.2044,0.009954,0.01319]

n = []

for i in range(0, len(tau_reg)):
    n = np.append(n, tau_reg[i] / (sigma_v * lcm[i]))
    
print('N = ', n)


# In[94]:


# MASS OF BOXES

mass = 1.67*10**-24
mass_correct = mass/0.7381
vol = lcm**3
print('Vol = ', vol)

# Mass of each box in grams
# n solved for in cell above
massgrams = n*vol*mass_correct

# Mass of each box in solar masses
sun = 1.988*10**33 # Mass of sun
solarmass = massgrams/sun


# In[95]:


# MASS OF WHOLE JET

# Volume of whole jet
npix = 762
width = ((arcpix * np.sqrt(npix)) / conv) * dist * (3.085677*10**18)   # [cm]
vol = width**3   # [cm^3]
print('VOLUME = ', vol, 'cm^3')

# Avg. density in sections 3A, 3F, 3G
avgden = (n[2] + n[3] + n[4]) / 3
print('AVG. DENS. = ', avgden)

massgrams = (avgden * vol * mass_correct) / sun
print('MASS IN GRAMS = ', massgrams)


# In[96]:


############# START OF ORIGINAL CODE WRITTEN BY JALYN KRAUSE ##############
############     EDITED DATE: 02/04/2020 - 02/06/2020        #############

# THE FOLLOWING CODE PLOTS TRANSMISSION (eqn. 5 in Lauren's Thesis) VS. BAND (V, I, J)


# In[97]:


# PLOTS FOR 1ST CHUNK OF DATA

# ARRAYS HOLD TRANSMISSION VALUES DERIVED FROM LAUREN'S PLOTS ABOVE

i1_v = [0.3221, 0.3345, 0.2844, 0.3298, 0.3554, 0.2876, 0.3708, 0.2942, 0.3089, 0.3433]
i1_i = [0.5603, 0.6047, 0.4897, 0.5459, 0.6287, 0.5270, 0.5826, 0.5520, 0.5572, 0.6138]
i1_j = [0.7325, 0.7639, 0.7007, 0.7358, 0.7444, 0.6948, 0.7413, 0.7496, 0.7607, 0.7985]

# EACH WAVELENGTH BAND (V = 555nm, I = 814nm, J = 1250nm)

v = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555]
i = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814]
j = [1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250]


# PLOTS FOR 1ST CHUNK OF DATA

fig = plt.figure(figsize = (11,8))

# Plots Lauren's data calculated from plots above
#plt.scatter(v, i1_v, edgecolor='black', label="V = [0.3221, 0.3345, 0.2844, 0.3298, 0.3554, 0.2876, 0.3708, 0.2942, 0.3089, 0.3433]")
#plt.scatter(i, i1_i, s=9.3, edgecolor='black', zorder=10,label="I = [0.5603, 0.6047, 0.4897, 0.5459, 0.6287, 0.5270, 0.5826, 0.5520, 0.5572, 0.6138]")
#plt.scatter(j, i1_j, edgecolor='black', label="J = [0.7325, 0.7639, 0.7007, 0.7358, 0.7444, 0.6948, 0.7413, 0.7496, 0.7607, 0.7985]")

# Normalized v = 0.3, Connects lines to see crossovers in transmission
plt.plot([v, i, j], [[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3], i1_i, i1_j], marker='o',  linestyle='dashed')
plt.title('Transmitted Light vs. Wavelength (Region 1)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.xlim((0.5,1.3))
#plt.ylim((0.,1.))
plt.legend(loc=4)
#fig.show()  # (USE TO DISPLY PLOTS)

#plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_1/region1_all.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME

#plt.show() # displays plot(s) on screen (NON JUPYTER)

print('V = ', i1_v, 'NORMALIZED TO 0.3')
print('I = ', i1_i)
print('J = ', i1_j)

#fig.show() #(USE TO DISPLAY PLOTS)

# PLOTS FOR 2ND CHUNK OF DATA

# ARRAYS HOLD TRANSMISSION VALUES DERIVED FROM LAUREN'S PLOTS ABOVE

i2_v = [0.3795, 0.3865, 0.3919, 0.3651, 0.3203, 0.3374, 0.2898, 0.3258, 0.3467]
i2_i = [0.7132, 0.7223, 0.6805, 0.6953, 0.6464, 0.6425, 0.6391, 0.6732, 0.6761]
i2_j = [0.8054, 0.8373, 0.8188, 0.8195, 0.7896, 0.8144, 0.8191, 0.8092, 0.7924]

# EACH WAVELENGTH BAND (V = 555nm, I = 814nm, J = 1250nm)

v = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555]
i = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814]
j = [1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250]

# PLOTS FIG

fig = plt.figure(figsize = (11,8))

# Plots data from Lauren's plots above
#plt.scatter(v, i2_v, edgecolor='black', label="V = [0.3795, 0.3865, 0.3919, 0.3651, 0.3203, 0.3374, 0.2898, 0.3258, 0.3467]")
#plt.scatter(i, i2_i, edgecolor='black', label="I = [0.7132, 0.7223, 0.6805, 0.6953, 0.6464, 0.6425, 0.6391, 0.6732, 0.6761]")
#plt.scatter(j, i2_j, edgecolor='black', label="J = [0.8054, 0.8373, 0.8188, 0.8195, 0.7896, 0.8144, 0.8191, 0.8092, 0.7924]")

#plt.legend(loc = 4)

# Normalized v = 0.3, Connects lines to see crossovers in transmission
plt.plot([v, i, j], [[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3], i2_i, i2_j], marker='o',  linestyle='dashed')


# In[ ]:


plt.title('Transmitted Light vs. Wavelength (Region 2)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.ylim((0.,1.))
#fig.show()  # (USE TO DISPLY PLOTS) 
#plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_2/region2_all.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME

print('V = ', i2_v, 'NORMALIZED TO 0.3')
print('I = ', i2_i)
print('J = ', i2_j)


# In[ ]:


# PLOTS FOR 3ND CHUNK OF DATA

i3_v = [0.2894, 0.3244, 0.3487, 0.4144, 0.4668, 0.3226, 0.2603, 0.2895, 0.2550, 0.3434]
i3_i = [0.4487, 0.5853, 0.6992, 0.8281, 0.8084, 0.4899, 0.4919, 0.5872, 0.4678, 0.5532]
i3_j = [0.6487, 0.7449, 0.8759, 0.9313, 0.8894, 0.6729, 0.6906, 0.7305, 0.7501, 0.7528]

# EACH WAVELENGTH BAND (V = 555nm, I = 814nm, J = 1250nm)

v = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555]
i = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814]
j = [1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250]

# PLOTTING FIG

fig = plt.figure(figsize = (11,8))

# Plots numbers from lauren's plots above
#plt.scatter(v, i3_v, edgecolor='black', label="V = [0.2894, 0.3244, 0.3487, 0.4144, 0.4668, 0.3226, 0.2603, 0.2895, 0.2550, 0.3434]")
#plt.scatter(i, i3_i, edgecolor='black', label="I = [0.4487, 0.5853, 0.6992, 0.8281, 0.8084, 0.4899, 0.4919, 0.5872, 0.4678, 0.5532]")
#plt.scatter(j, i3_j, edgecolor='black', label="J = [0.6487, 0.7449, 0.8759, 0.9313, 0.8894, 0.6729, 0.6906, 0.7305, 0.7501, 0.7528]")
#plt.legend(loc = 4)


# Normalized v = 0.3, Connects lines to see crossovers in transmission
plt.plot([v, i, j], [[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3], i3_i, i3_j], marker='o',  linestyle='dashed')

plt.title('Transmitted Light vs. Wavelength (Region 3)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.ylim((0.,1.))
#fig.show()  # (USE TO DISPLY PLOTS)
#plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_3/region3_all.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME


print('V = ', i3_v, 'NORMALIZED 0.3')
print('I = ', i3_i)
print('J = ', i3_j)


# In[ ]:


# PLOTS ALL DATA FROM THE TRANSMISSION DATA FOR EACH REGION ON THE SAME DIAGRAM

fig = plt.figure(figsize = (11,8))

# XDATA 
v_ten = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555]  # ten elements in array
v_nine = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555]  # nine elements in array
i_ten = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814]
i_nine = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814]
j_ten = [1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250]
j_nine = [1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250]


r_one = plt.plot((v_ten,i_ten,j_ten),(i1_v,i1_i,i1_j), marker='o',ms=4, lw=1.25, color='b',  linestyle='dashed')
r_two = plt.plot((v_nine,i_nine,j_nine),(i2_v,i2_i,i2_j), marker='o',ms=4, lw=1.25, color='r', linestyle='dashdot')
r_three = plt.plot((v_ten,i_ten,j_ten),(i3_v,i3_i,i3_j), marker='o',ms=4, lw=1.25, color='g',  linestyle='dotted')

line1, = plt.plot([1, 2, 3], label="REGION 1 -- DASHED BLUE LINE", color='b', linestyle='dashed')
line2, = plt.plot([3, 2, 1], label="REGION 2 -- DASHED DOT RED LINE", color='r', linestyle='dashdot')
line3, = plt.plot([2, 2, 2], label="REGION 3 -- DOTTED GREEN LINE", color='g', linestyle='dotted')

plt.title('Transmitted Light vs. Wavelength (All Regions)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.xlim((0.5,1.3))
plt.ylim((0.,1.))
plt.legend(loc=4)
#fig.show() # (USE TO DISPLAY PLOTS)
#plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/all_trans_data.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME

print('V_one = ', i1_v)
print('I_one = ', i1_i)
print('J_one = ', i1_j)

print('V_two = ', i2_v)
print('I_two = ', i2_i)
print('J_two = ', i2_j)

print('V_three = ', i3_v)
print('I_three = ', i3_i)
print('J_three = ', i3_j)

#plt.close("all") # CLOSES ALL PLOTS FROM THE PREVIOUS SECTIONS

# In[ ]:

#########################################################################################
#########################################################################################


# Edited by Jalyn Krause on 20 May 2020 - PRESENT

# USES UNIFORM COLOR MODEL (APPROACH 2)


# In[ ]:


# We assume that the galaxy has a uniform color.
# Thus we are saying, I0(lambda1, x, y)/I0(lambda2, x, y) ~ constant = Iratio(lambda1,lambda2).
# We use our assumption that interstellar obscuration increases monotonically as wavelength decreases.
# In the simplest case, as an example, if a dust cloud is present with optical depth tau(V), I (lambda) = I0(lambda)exp[-k(lambda)*tauV].
# If we take the ratio of 2 images that are normalized such that the ratio for unobscured stars is set equal to 1, 
#      i.e,. Iratio0(lambda1,lambda2)=1 then in an obscured location we have TrColor = I(ratio, lambda1, lambda2) = Iratio0 * exp[-(k(lambda1)-k(lambda2)*tauV]
# Since k(lambda) is known for our extinction curve, we can derive tauV for the cloud.
# Our data provide TrColor for HST V/I and V/J, and I/J. 
# We also have HST U/V. It's low signal-to-noise but allows us to find lower levels of obscuration.
# HST is diffraction limited, so the angular resolution scales as ~1.22 lambda/(aperture diameter) ~0.04 (lambda/555 nm) arcsec for lambda >500 nm. 
# Thus each image has different resolution. 
# In addition, the IR F125W image was obtain with a different detector and has larger pixels than the optical CCDs in WFC3.
# Ralf resampled the images to a common pixel scale using astrodrizzle. 
# He then smoothed the shorter wavelength images so that their resolution would match that of the lowest angular resolution F125W image.
# The resampled and psf-corrected TrColor images are the zw049*trans_ratio*psf.fits in UW Box in a new directory, ColorMaps.
# The zw049*_norm files contain ratios of images where the unobscured stellar background ratio has been set to 1.0+/-~0.05. These ten give us multiple estimates of TrColor.


# In[ ]:

''' 
# PLOTS TRANSMISSION RATIOS FOR ALL REGIONS


# REGION 1, RATIO TO I
ri_vi_rat = [v / x for v,x in zip(i1_v,i1_i)]
ri_ii_rat = [i / x for i,x in zip(i1_i,i1_i)]
ri_ji_rat = [j / x for j,x in zip(i1_j,i1_i)]

# REGION 1, RATIO TO J
ri_vj_rat = [v / y for v,y in zip(i1_v,i1_j)]
ri_ij_rat = [i / y for i,y in zip(i1_i,i1_j)]
ri_jj_rat = [j / y for j,y in zip(i1_j,i1_j)]

# REGION 1, RATIO TO V
ri_vv_rat = [v / z for v,z in zip(i1_v,i1_v)]
ri_iv_rat = [i / z for i,z in zip(i1_i,i1_v)]
ri_jv_rat = [j / z for j,z in zip(i1_j,i1_v)]

# PLOTTING FIG FOR REGION 1
r_one_i = plt.plot((v_ten,i_ten,j_ten),(ri_vi_rat,ri_ii_rat,ri_ji_rat), marker='o',ms=4, lw=1.25, color='b',  linestyle='dashed')
r_one_j = plt.plot((v_ten,i_ten,j_ten),(ri_vj_rat,ri_ij_rat,ri_jj_rat), marker='o',ms=4, lw=1.25, color='r', linestyle='dashdot')
r_one_v = plt.plot((v_ten,i_ten,j_ten),(ri_vv_rat,ri_iv_rat,ri_jv_rat), marker='o',ms=4, lw=1.25, color='g',  linestyle='dotted')

plt.title('Transmitted Light vs. Wavelength (Region 1)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.xlim((0.5,1.3))
plt.legend(loc=4)
#plt.show()  # (USE TO DISPLY PLOTS)
'''

############################################################################################################################################################
# BEGINING OF CODE FOR PLOTTING DATA FOR INTERSECTION POINT BETWEEN V,I, AND J BANDS FOR ADDITIONAL REGIONS (POLAR RING) DEFINED BY JALYN KRAUSE IN AUG. 2020


# Sections 4A-J of regions defined by Jalyn K.
# Diagram of sections shown in data excel file 

# Values taken from ds9 data for 'Nick's Pillar' regions defined by Lauren L.
v_4 = [0.6899, 0.7317, 0.5302, 0.6104, 0.6210, 0.6353, 0.5991, 0.6573, 0.6046, 0.5473]
i_4 = [0.8204, 0.8330, 0.6516, 0.6419, 0.7448, 0.6755, 0.6558, 0.8259, 0.6504, 0.6141]
j_4 = [0.8327, 0.8595, 0.7786, 0.7245, 0.8259, 0.7205, 0.6847, 0.8128, 0.7051, 0.7593]

# Developed from equn (5) of Lauren's senior thesis
# V = f + (1-f)*e^(-tau_V)
# ==> tau_v = log( (V -f) / 1-f)
# f = fractional distance into the cloud where the dust is located
# Left side of equation is a ratio of the total light the observer sees in each band, in a fractional form compared to the total light emitted
# *Recall all tau values are in terms of tau_v*

def tau(f,V):
    xxx = (V-f) / (1-f)
    yyy = -1*np.log( xxx )
    yyy[xxx <= 0] = 5
    return yyy

f_input = np.linspace(0,1,101)

filenames = ['4A', '4B',  '4C', '4D', '4E', '4F', '4G', '4H', '4I', '4J']
    
for i in range(0, len(filenames)):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(f_input, tau(f_input, v_4[i])/3.1*3.1, label="V="+str(v_4[i]))
    ax.plot(f_input, tau(f_input, i_4[i])/1.87*3.1, label="I="+str(i_4[i]))
    ax.plot(f_input, tau(f_input, j_4[i])/0.87*3.1, label="J="+str(j_4[i]))
    ax.set_ylim((0.,3))
    ax.set_xlim((0,1))
    ax.set_title(r'$\tau_v$ vs. Fractional Location in Galaxy (f) (REGION 4)')
    ax.set_xlabel('f')
    ax.set_ylabel(r'$\tau_v$')
    ax.legend(loc = 'best')
    #print('V = ', v_4[i], 'I = ', i_4[i], 'J = ', j_4[i])
        
    #plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_4/'+ filenames[i] + '.png') #SAVES EACH FIG TO FOLDER LABELED BY REGION NAME
    
    
# Sections 5A-J of regions defined by Jalyn K.
# Diagram of sections shown in data excel file 

v_5 = [0.817346, 0.851550, 0.731123, 0.816211, 0.764386, 0.905439, 0.821342, 0.822922, 0.877841, 0.934868]
i_5 = [0.762223, 0.741940, 0.738229, 0.793757, 0.762880, 0.797480, 0.810254, 0.829480, 0.814794, 0.859911]
j_5 = [0.830558, 0.864141, 0.850765, 0.895060, 0.875499, 0.911799, 0.882971, 0.914812, 0.903250, 0.912751]

# Developed from equn (5) of Lauren's senior thesis
# V = f + (1-f)*e^(-tau_V)
# ==> tau_v = log( (V -f) / 1-f)
# f = fractional distance into the cloud where the dust is located
# Left side of equation is a ratio of the total light the observer sees in each band, in a fractional form compared to the total light emitted
# *Recall all tau values are in terms of tau_v*

def tau(f,V):
    xxx = (V-f) / (1-f)
    yyy = -1*np.log( xxx )
    yyy[xxx <= 0] = 5
    return yyy

f_input = np.linspace(0,1,101)

filenames = ['5A', '5B',  '5C', '5D', '5E', '5F', '5G', '5H', '5I', '5J']
    
for j in range(0, len(filenames)):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(f_input, tau(f_input, v_5[j])/3.1*3.1, label="V="+str(v_5[j]))
    ax.plot(f_input, tau(f_input, i_5[j])/1.87*3.1, label="I="+str(i_5[j]))
    ax.plot(f_input, tau(f_input, j_5[j])/0.87*3.1, label="J="+str(j_5[j]))
    ax.set_ylim((0.,3))
    ax.set_xlim((0,1))
    ax.set_title(r'$\tau_v$ vs. Fractional Location in Galaxy (f) (REGION 5)')
    ax.set_xlabel('f')
    ax.set_ylabel(r'$\tau_v$')
    ax.legend(loc = 'best')
    #print('V = ', v_5[j], 'I = ', i_5[j], 'J = ', j_5[j])
    
    #plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_5/'+ filenames[j] + '.png') #SAVES EACH FIG TO FOLDER LABELED BY REGION NAME

# Sections 6A-J of regions defined by Jalyn K.
# Diagram of sections shown in data excel file 

v_6 = [0.866715, 0.918022, 0.663775, 0.550166, 0.753511, 0.640885, 0.542817, 0.571696, 0.747268, 0.710120]
i_6 = [0.925089, 0.879423, 0.760106, 0.614777, 0.730369, 0.706440, 0.622697, 0.623209, 0.690437, 0.639051]
j_6 = [0.909783, 0.857721, 0.914602, 0.842134, 0.815482, 0.926424, 0.881680, 0.845964, 0.837974, 0.889399]

# Developed from equn (5) of Lauren's senior thesis
# V = f + (1-f)*e^(-tau_V)
# ==> tau_v = log( (V -f) / 1-f)
# f = fractional distance into the cloud where the dust is located
# Left side of equation is a ratio of the total light the observer sees in each band, in a fractional form compared to the total light emitted
# *Recall all tau values are in terms of tau_v*

def tau(f,V):
    xxx = (V-f) / (1-f)
    yyy = -1*np.log( xxx )
    yyy[xxx <= 0] = 5
    return yyy

f_input = np.linspace(0,1,101)

filenames = ['6A', '6B',  '6C', '6D', '6E', '6F', '6G', '6H', '6I', '6J']
    
for k in range(0, len(filenames)):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(f_input, tau(f_input, v_6[k])/3.1*3.1, label="V="+str(v_6[k]))
    ax.plot(f_input, tau(f_input, i_6[k])/1.87*3.1, label="I="+str(i_6[k]))
    ax.plot(f_input, tau(f_input, j_6[k])/0.87*3.1, label="J="+str(j_6[k]))
    ax.set_ylim((0.,3))
    ax.set_xlim((0,1))
    ax.set_title(r'$\tau_v$ vs. Fractional Location in Galaxy (f) (REGION 6)')
    ax.set_xlabel('f')
    ax.set_ylabel(r'$\tau_v$')
    ax.legend(loc = 'best')
    #print('V = ', v_6[k], 'I = ', i_6[k], 'J = ', j_6[k])
    
    #plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_6/'+ filenames[k] + '.png') #SAVES EACH FIG TO FOLDER LABELED BY REGION NAME

# Sections 7A-K of regions defined by Jalyn K.
# Diagram of sections shown in data excel file 

v_7 = [0.831520, 0.811042, 0.747706, 0.740157, 0.680097, 0.729802, 0.688435, 0.641984, 0.702126, 0.811180, 0.819490]
i_7 = [0.788054, 0.743663, 0.683360, 0.679449, 0.679134, 0.703320, 0.641194, 0.605023, 0.653657, 0.680700, 0.650440]
j_7 = [0.832171, 0.815342, 0.787296, 0.831368, 0.798368, 0.814670, 0.802350, 0.869040, 0.837509, 0.851313, 0.845949]

# Developed from equn (5) of Lauren's senior thesis
# V = f + (1-f)*e^(-tau_V)
# ==> tau_v = log( (V -f) / 1-f)
# f = fractional distance into the cloud where the dust is located
# Left side of equation is a ratio of the total light the observer sees in each band, in a fractional form compared to the total light emitted
# *Recall all tau values are in terms of tau_v*

def tau(f,V):
    xxx = (V-f) / (1-f)
    yyy = -1*np.log( xxx )
    yyy[xxx <= 0] = 5
    return yyy

f_input = np.linspace(0,1,101)

filenames = ['7A', '7B',  '7C', '7D', '7E', '7F', '7G', '7H', '7I', '7J', '7K']
    
for k in range(0, len(filenames)):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(f_input, tau(f_input, v_7[k])/3.1*3.1, label="V="+str(v_7[k]))
    ax.plot(f_input, tau(f_input, i_7[k])/1.87*3.1, label="I="+str(i_7[k]))
    ax.plot(f_input, tau(f_input, j_7[k])/0.87*3.1, label="J="+str(j_7[k]))
    ax.set_ylim((0.,3))
    ax.set_xlim((0,1))
    ax.set_title(r'$\tau_v$ vs. Fractional Location in Galaxy (f) (REGION 7)')
    ax.set_xlabel('f')
    ax.set_ylabel(r'$\tau_v$')
    ax.legend(loc = 'best')
    #print('V = ', v_7[k], 'I = ', i_7[k], 'J = ', j_7[k])
    
    #plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_7/'+ filenames[k] + '.png') #SAVES EACH FIG TO FOLDER LABELED BY REGION NAME



# PLOTS ALL DATA FOR EACH REGION ON ONE PLOT

# EACH WAVELENGTH BAND (V = 555nm, I = 814nm, J = 1250nm)

v = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555]
i = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814]
j = [1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250]

# PLOTS FOR 4th CHUNK OF DATA

fig = plt.figure(figsize = (11,8))

# Normalized v = 0.3, Connects lines to see crossovers in transmission
plt.plot([v, i, j], [[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3], i_4, j_4], marker='o',  linestyle='dashed')

plt.title('Transmitted Light vs. Wavelength (Region 4)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.ylim((0.,1.))
#fig.show()  # (USE TO DISPLY PLOTS)
#plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_4/region4_all.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME


print('V = ', v_4, 'NORMALIZED 0.3')
print('I = ', i_4)
print('J = ', j_4)

# PLOTS FOR 5th CHUNK OF DATA

fig = plt.figure(figsize = (11,8))

# Normalized v = 0.3, Connects lines to see crossovers in transmission
plt.plot([v, i, j], [[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3], i_5, j_5], marker='o',  linestyle='dashed')
plt.title('Transmitted Light vs. Wavelength (Region 5)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.ylim((0.,1.))
plt.legend(loc=4)
#fig.show()  # (USE TO DISPLY PLOTS)
#plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_5/region5_all.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME

print('V = ', v_5, 'NORMALIZED TO 0.3')
print('I = ', i_5)
print('J = ', j_5)

# PLOTS FOR 6th CHUNK OF DATA
fig = plt.figure(figsize = (11,8))

# Normalized v = 0.3, Connects lines to see crossovers in transmission
plt.plot([v, i, j], [[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3], i_6, j_6], marker='o',  linestyle='dashed')
plt.title('Transmitted Light vs. Wavelength (Region 6)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.ylim((0.,1.))
#fig.show()  # (USE TO DISPLY PLOTS) 
#plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_6/region6_all.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME

print('V = ', v_6, 'NORMALIZED TO 0.3')
print('I = ', i_6)
print('J = ', j_6)


# PLOTS FOR 7th CHUNK OF DATA

# EACH WAVELENGTH BAND (V = 555nm, I = 814nm, J = 1250nm)

v = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555, 0.555]
i = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814, 0.824]
j = [1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250, 1.250]


fig = plt.figure(figsize = (11,8))

# Normalized v = 0.3, Connects lines to see crossovers in transmission
plt.plot([v, i, j], [[0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3], i_7, j_7], marker='o',  linestyle='dashed')

plt.title('Transmitted Light vs. Wavelength (Region 7)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.ylim((0.,1.))
#fig.show()  # (USE TO DISPLY PLOTS)
#plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/region_7/region7_all.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME


print('V = ', v_7, 'NORMALIZED 0.3')
print('I = ', i_7)
print('J = ', j_7)


# In[ ]:


# PLOTS ALL DATA FROM THE TRANSMISSION DATA FOR EACH REGION ON THE SAME DIAGRAM

fig = plt.figure(figsize = (11,8))

# XDATA 
v_ten = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555]  # ten elements in array
v_el = [0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555,0.555]  # eleven elements in array
i_ten = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814]
i_el = [0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814,0.814]
j_ten = [1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250]
j_el = [1.250, 1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250,1.250]


r_four = plt.plot((v_ten,i_ten,j_ten),(v_4,i_4,j_4), marker='o',ms=4, lw=1.25, color='b',  linestyle='dashed')
r_five = plt.plot((v_ten,i_ten,j_ten),(v_5,i_5,j_5), marker='o',ms=4, lw=1.25, color='r', linestyle='dashdot')
r_six = plt.plot((v_ten,i_ten,j_ten),(v_6,i_6,j_6), marker='o',ms=4, lw=1.25, color='g',  linestyle='dotted')
r_seven = plt.plot((v_el,i_el,j_el),(v_7,i_7,j_7), marker='o',ms=4, lw=1.25, color='m',  linestyle='dotted')

line1, = plt.plot([1, 2, 3], label="REGION 4 -- DASHED BLUE LINE", color='b', linestyle='dashed')
line2, = plt.plot([3, 2, 1], label="REGION 5 -- DASHED DOT RED LINE", color='r', linestyle='dashdot')
line3, = plt.plot([2, 2, 2], label="REGION 6 -- DOTTED GREEN LINE", color='g', linestyle='dotted')
line4, = plt.plot([2, 2, 2], label="REGION 7 -- SOLID MAGENTA LINE", color='m', linestyle='solid')

plt.title('Transmitted Light vs. Wavelength (All Regions 4-7)')
plt.xlabel('Wavelength (Microns) (V, I, J bands)')
plt.ylabel('Trans.')
plt.xlim((0.5,1.3))
plt.ylim((0.,1.))
plt.legend(loc=4)
#fig.show() # (USE TO DISPLAY PLOTS)
plt.savefig('/home/krause/Documents/gallagher/reu_2020/lauren_plots/all_trans_data_4-7.png') #SAVES FIG TO FOLDER LABELED BY REGION NAME

'''
print('V_one = ', i1_v)
print('I_one = ', i1_i)
print('J_one = ', i1_j)

print('V_two = ', i2_v)
print('I_two = ', i2_i)
print('J_two = ', i2_j)

print('V_three = ', i3_v)
print('I_three = ', i3_i)
print('J_three = ', i3_j)
'''

