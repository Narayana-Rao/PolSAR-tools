
""" 
This script generates H-Theta_FP clustering output for Full-polarimetric data.

Implemented from Dey, S., Bhattacharya, A., Ratha, D., Mandal, D., McNairn, H., Lopez-Sanchez, J.M. and Rao, Y.S., 2020. 
"Novel clustering schemes for full and compact polarimetric SAR data: An application for rice phenology characterization." 
ISPRS Journal of Photogrammetry and Remote Sensing, 169, pp.135-151.

"""
from osgeo import gdal
import numpy as np
from tkinter import filedialog
from tkinter import *
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap


def read_bin(file):
    ds = gdal.Open(file)
    gt=ds.GetGeoTransform()
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()
    return arr,gt

def write_bin(file,wdata,refData):
            
    ds = gdal.Open(refData)
    [cols, rows] = wdata.shape
            
    driver = gdal.GetDriverByName("ENVI")
    outdata = driver.Create(file, rows, cols, 1, gdal.GDT_Float32)
    outdata.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
    outdata.SetProjection(ds.GetProjection())##sets same projection as input
                
    outdata.SetDescription(file)
    outdata.GetRasterBand(1).WriteArray(wdata)
    outdata.FlushCache() 

root = Tk()
root.withdraw()
folder = filedialog.askdirectory()

[theta_dat, gt1] = read_bin(folder + '/Theta_FP.bin')
[ent_dat, gt2] = read_bin(folder + '/H_FP.bin')

ent_dat = 1 - ent_dat

cluster_img = np.ones(np.shape(theta_dat))

cluster_img[(ent_dat > 0.5) & (ent_dat <= 1) & (theta_dat >= -90) & (theta_dat < -10)] = 1
cluster_img[(ent_dat > 0.5) & (ent_dat <= 1) & (theta_dat >= -10) & (theta_dat < 0)] = 4
cluster_img[(ent_dat > 0.5) & (ent_dat <= 1) & (theta_dat >= 0) & (theta_dat < 20)] = 7
cluster_img[(ent_dat > 0.5) & (ent_dat <= 1) & (theta_dat >= 20) & (theta_dat <= 90)] = 10

cluster_img[(ent_dat > 0.3) & (ent_dat <= 0.5) & (theta_dat >= -90) & (theta_dat < -10)] = 2
cluster_img[(ent_dat > 0.3) & (ent_dat <= 0.5) & (theta_dat >= -10) & (theta_dat < 0)] = 5
cluster_img[(ent_dat > 0.3) & (ent_dat <= 0.5) & (theta_dat >= 0) & (theta_dat < 20)] = 8
cluster_img[(ent_dat > 0.3) & (ent_dat <= 0.5) & (theta_dat >= 20) & (theta_dat <= 90)] = 11

cluster_img[(ent_dat >= 0) & (ent_dat <= 0.3) & (theta_dat >= -90) & (theta_dat < -10)] = 3
cluster_img[(ent_dat >= 0) & (ent_dat <= 0.3) & (theta_dat >= -10) & (theta_dat < 0)] = 6
cluster_img[(ent_dat >= 0) & (ent_dat <= 0.3) & (theta_dat >= 0) & (theta_dat < 20)] = 9
cluster_img[(ent_dat >= 0) & (ent_dat <= 0.3) & (theta_dat >= 20) & (theta_dat <= 90)] = 12

cluster_img = cluster_img.astype(np.float32)
cluster_img[np.isnan(theta_dat) == True] = np.nan


mymap = [
    (1,0.69,0.69),
    (0.97, 0.24, 0.24),
    (0.77, 0.19, 0.19),
    (0.51,1,0.84),
    (0.01, 1, 0.67),
    (0.50, 0.99, 0.51),
    (0.45, 0.89, 0.46),
    (0.35, 0.69, 0.58),
    (0.11,0.44,0.23),
    (0.51,0.58,1),
    (0.2,0.32,0.97),
    (0,0.15,0.96)
   ]

mymap = LinearSegmentedColormap.from_list("mymap", mymap, N=12)


num_list = []

low_val = 1
high_val = 12

step = (high_val - low_val)/24 # total 12 zones
cur_val = low_val + step
num_list.append(cur_val)

for i in np.arange(low_val+1, high_val+1):
    step1 = 2*step 
    cur_val = cur_val + step1
    num_list.append(cur_val)

im_ratio = cluster_img.shape[0]/cluster_img.shape[1]

fig = plt.figure(figsize=(10,8))
plt.imshow(cluster_img, cmap = mymap)
cbar = mpl.pyplot.colorbar(fraction=0.047*im_ratio)

cbar.set_ticks(num_list)
cbar.set_ticklabels(["Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z9", "Z10", "Z11", "Z12"])
for t in cbar.ax.get_yticklabels():
     t.set_fontsize(20)

save_path = folder + '/cluster_FP.png'
fig.savefig(save_path,dpi=300, bbox_inches='tight')
#%%
ofilegrvi = folder +'/cluster_FP.bin'
infile = folder + '/Theta_FP.bin'
write_bin(ofilegrvi,cluster_img,infile)

#%% Data extraction

from math import floor
from osgeo import ogr

# shp_filename = 'test_vector.shp'


# ds=ogr.Open(shp_filename)
# lyr=ds.GetLayer()

# points = []

# for feat in lyr:
#     geom = feat.GetGeometryRef()
#     mx,my=geom.GetX(), geom.GetY()  #coord in map units

#     #Convert from map to pixel coordinates.
#     #Only works for geotransforms with no rotation.
#     px = floor((mx - gt1[0]) / gt1[1]) #x pixel
#     py = floor((my - gt1[3]) / gt1[5]) #y pixel
#     points.append((py,px))

# pt = np.array(points)


# # get the data
# theta_points = theta_dat[pt[:,0],pt[:,1]]
# ent_points = ent_dat[pt[:,0],pt[:,1]]  
# cluster_img_points = cluster_img[pt[:,0],pt[:,1]]

theta_points = theta_dat.flatten()
ent_points = ent_dat.flatten()  

#%% Plotting-----
plt.rcParams.update({'font.size': 8})
def flag_df(df):

    if (df['H'] >0.5) and  (df['Theta_DP'] >=-90) and (df['Theta_DP'] < -10):
        return '1'
    elif (df['H'] >=0.3) and (df['H'] <0.5) and (df['Theta_DP'] >=-90) and (df['Theta_DP'] < -10):
        return '2'
    elif (df['H'] >=0) and (df['H'] <0.3) and  (df['Theta_DP'] >=-90) and (df['Theta_DP'] < -10):
        return '3'
    
    elif (df['H'] >=0) and (df['H'] <0.3) and (df['Theta_DP'] <0) and (df['Theta_DP'] >=-10):
        return '6'
    elif (df['H'] >=0.3) and (df['H'] <0.5) and (df['Theta_DP'] <0) and (df['Theta_DP'] >= -10):
        return '5'
    elif (df['H'] >=0.5) and (df['H'] <1) and (df['Theta_DP'] <0) and (df['Theta_DP'] >= -10):
        return '4'
    
    elif (df['H'] >=0) and (df['H'] <0.3) and (df['Theta_DP'] >=0) and (df['Theta_DP'] <20):
        return '9'
    elif (df['H'] >=0.3) and (df['H'] <0.5) and (df['Theta_DP'] >=0) and (df['Theta_DP'] <20):
        return '8'
    elif (df['H'] >=0.5) and (df['H'] <1) and (df['Theta_DP'] >=0) and (df['Theta_DP'] <20):
        return '7'    

    elif (df['H'] >=0) and (df['H'] <0.3) and (df['Theta_DP'] >=20) and (df['Theta_DP'] <=90):
        return '12'
    elif (df['H'] >=0.3) and (df['H'] <0.5) and (df['Theta_DP'] >=20) and (df['Theta_DP'] <=90):
        return '11'
    elif (df['H'] >=0.5) and (df['H'] <1) and (df['Theta_DP'] >=20) and (df['Theta_DP'] <=90):
        return '10'    

    
    else:
        return 0


import pandas as pd

plot_df = pd.DataFrame()
plot_df['H'] =ent_points
plot_df['Theta_DP'] = theta_points

plot_df['C'] = plot_df.apply(flag_df, axis = 1)
lw=0.3
fig = plt.figure(figsize=(3,2.3),dpi=300)
ax = plt.subplot(111, polar=True)
# set zero north
ax.set_theta_zero_location('N')
ax.set_theta_direction('clockwise')


# Theoretical boundary
# m=0
c1l=[]
for m in np.arange(0,1.005,0.005):
    T31 = np.array([[1,0,0],[0,m,0],[0,0,m]])   
    mfp = np.sqrt(1-((27*np.linalg.det(T31))/(np.trace(T31))**3))
    tfp1 = 2*np.arctan(mfp*np.trace(T31)*(1-m-m)/(1*(m+m)+np.trace(T31)**2*mfp**2))*180/np.pi
    
    eigenValues, eigenVectors = np.linalg.eig(T31)
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    
    p1 = eigenValues[0]/eigenValues.sum()
    p2 = eigenValues[1]/eigenValues.sum()
    p3 = eigenValues[2]/eigenValues.sum()
            
    h1=-(p1*np.log(p1) / np.log(3)+ p2*np.log(p2) / np.log(3) +p3*np.log(p3) / np.log(3))
    
    if tfp1==90:
        h1=0
        
    c1l.append([tfp1,1-h1,mfp])

# curve-I
ax.plot((np.array(c1l)[:,0]*np.pi/180),np.array(c1l)[:,1],'k-',linewidth = lw)
ax.fill_between((np.array(c1l)[:,0]*np.pi/180),np.array(c1l)[:,1],color='gray', alpha=0.6,edgecolor="b", linewidth=0.0)
c2l=[[-90,1]]
c2l=[]

for m in np.arange(0.5,1.005,0.005):
    T31 = np.array([[2*m-1,0,0],[0,1,0],[0,0,1]])    
    mfp = np.sqrt(1-((27*np.linalg.det(T31))/(np.trace(T31))**3))
    tfp1 = 2*np.arctan(mfp*np.trace(T31)*(2*m-1-1-1)/((2*m-1)*(1+1)+np.trace(T31)**2*mfp**2))*180/np.pi
    
    eigenValues, eigenVectors = np.linalg.eig(T31)
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    
    p1 = eigenValues[0]/eigenValues.sum()
    p2 = eigenValues[1]/eigenValues.sum()
    p3 = eigenValues[2]/eigenValues.sum()
            
    h1=-(p1*np.log(p1) / np.log(3)+ p2*np.log(p2) / np.log(3) +p3*np.log(p3) / np.log(3))
           
    c2l.append([tfp1,1-h1,mfp])



# curve-II
c2l = np.array(c2l)
c2l = c2l[~np.isnan(c2l).any(axis=1)]

ax.plot((np.array(c2l)[:,0]*np.pi/180),np.array(c2l)[:,1],'k-',linewidth = lw)
ax.fill_between((np.array(c2l)[:,0]*np.pi/180),np.array(c2l)[:,1],color='gray', alpha=0.6, edgecolor="b", linewidth=0.0)

#Curve-III
ax.plot(np.repeat(-90*np.pi/180,np.size(np.arange(0.34,1.1,.1))),np.arange(0.34,1.1,.1),'k-',linewidth = lw)

ax.plot((0,20*np.pi/180),(0,1),'k--',linewidth = lw)
ax.plot((0,-10*np.pi/180),(0,1),'k--',linewidth = lw)
ax.plot((0,0),(0,1),'k--',linewidth = lw)

theta = np.linspace(-np.pi,np.pi)
r = np.repeat(0.3, np.size(theta), axis=None)
ax.plot(theta, r, 'k--',linewidth = lw)
r = np.repeat(0.5, np.size(theta), axis=None)
ax.plot(theta, r, 'k--',linewidth = lw)


cmap = mpl.colors.ListedColormap([#
                                  
#'#000000',
'#ffc4b2',
'#ff544a',
'#c90305',

'#efffff',
'#03fdc5',
'#b2fb9a',

'#a4ff86',
'#82daa8',
'#2c964a',

'#b7bffe', #Z10
'#4e71ff', #Z11
'#012fff', #Z12
])
    
bounds = [ 1, 2, 3, 4,5,6,7,8,9,10,11,12,13]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

ax.scatter(plot_df['Theta_DP']*np.pi/180,plot_df['H'],label='Observed',c=np.array(plot_df['C'],dtype=int),
           cmap=cmap,vmin=1,vmax=12,norm=norm,alpha=1,s=0.5)

ax.text(0.44, 0.2, 
      '  'r'{:.1f}' '   {:.1f}'' {:.1f}''        {:.1f}'.format(0.0,0.3,0.5,1.0),
        transform=ax.transAxes,
        )

ax.text(0.7, 0.12, 
      '  'r'$\overline{H}$',
        transform=ax.transAxes,
        )

ax.text(0.45, 0.9, 
      '  'r'$\overline{\theta}_{FP}$',
        transform=ax.transAxes,
        )

ax.set_ylim(0,1)
ax.xaxis.set_tick_params(pad=0) 
ax.set_yticks(np.array([]))
ax.set_xticks(np.array([-90, -10, 0, 20, 90])/180*np.pi)
ax.set_thetalim(-1/2*np.pi, 1/2*np.pi)


for side in ax.spines.keys(): 
    ax.spines[side].set_linewidth(0.2)

plt.grid(False)
ax.annotate('',
            xy=(0.37,0.88), xycoords='axes fraction',
            xytext=(0.47,.91), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="angle3,angleA=0,angleB=-150",
                            color='k',
                            linewidth=0.3))
   
ax.annotate('',
            xy=(0.66,0.88), xycoords='axes fraction',
            xytext=(0.56,0.91), textcoords='axes fraction',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="angle3,angleA=0,angleB=-30",
                            color='k',
                            linewidth=0.3)) 


plt.tight_layout()

save_path = folder + '/h_theta_FP.png'
plt.savefig(save_path,dpi=400)