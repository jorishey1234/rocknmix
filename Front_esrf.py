#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import argparse
import datetime
import time
import sys

timestamp = datetime.datetime.fromtimestamp(time.time()).strftime("%Y-%m-%d %H:%M:%S")

print("="*(65))
print('MIXING FRONT ANALYSIS - '+ timestamp)
print("="*(65))

parser = argparse.ArgumentParser(description='Front analysis')
parser.add_argument('-S','--SAMPLE', default='./',help='Sample name')
parser.add_argument('-R','--RESOLUTION', default='LOW',help='Scan resolution(HIGH/LOW/Med)')
parser.add_argument('-T','--TYPE', default='COFLOW',help='Experiment type (COFLOW_Q4)')
parser.add_argument('-path',default='./',help='Path of result directory')
parser.add_argument('-hdf5name',default='I',help='Name of image in hdf5 file')
parser.add_argument('-date',default=timestamp,help='Time stamp')
parser.add_argument('-rot',default=0,help='Set rotation of the front (0 automatic)')

parser.add_argument('-reg', default=[],help='Location for image registration between dry and wet. Sould be outside of the front')
parser.add_argument('-c0', default=0.0,help='min intensity of the clear fluid')
parser.add_argument('-c1', default=1.0,help='max intensity of the brine')
parser.add_argument('-dryth', default=0.05,help='Threshold for solid segmentation in the dry scan')
parser.add_argument('-nstep', default=50,help='Number of computing planes in the downstream direction')
parser.add_argument('-nmin', default=0)
parser.add_argument('-nmax', default=-1)
parser.add_argument('-rad', default=950,help='Radius in px of the ROI')
parser.add_argument('-nbin', default=100,help='Number of bins for front statistics')

args = parser.parse_args()


#%%
# ================= Parameters =========================================
#fmt='tif' # image format
path=args.path

file_dry=path+'/'+args.SAMPLE+'/'+args.SAMPLE+'_'+'DRY'+'_'+args.RESOLUTION +'.hdf5' #Dry scan
file_wet=path+'/'+args.SAMPLE+'/'+args.SAMPLE+'_'+args.TYPE+'_'+args.RESOLUTION +'.hdf5' #W et scan

import h5py

print("Reading dry file at "+file_dry)
try:
	File_dry=h5py.File(file_dry,'r')
	print('Read ',File_dry[args.hdf5name].shape)
except:
	print("No file found")
	sys.exit()


print("Reading wet file at "+file_wet)

try:
	File_wet=h5py.File(file_wet,'r')
	print('Read ',File_wet[args.hdf5name].shape)
except:
	print("No file found")
	sys.exit()

# Image prefix, followed by 0001.tif, 0002.tif ...
#prefix_dry='dry-scan 24May'
#prefix_wet='TwoFluids_5_5mlmin_24th May'

#Parameters
nmin=args.nmin
nmax=args.nmax
if nmax==-1:# read all frames
	nmax=File_wet[args.hdf5name].shape[-1]
nstep=args.nstep



# Other parameters
p1=args.c0 # Approx concentration of clear fluid
p2=args.c1 # Approx concentration of brine
th=args.dryth # Concentration threshold for solid detection
rad=args.rad # Radius of core sample
reg=args.reg


dn=args.nbin # number of bins for study of front variance
ROTATE=2 # method for front rotation

# Image subset for registration
if len(reg)>0:
	sX=[reg[0]-100,reg[0]+100]
	sY=[reg[0]-100,reg[0]+100]

# ================= Parameters =========================================


import numpy as np
import matplotlib.pyplot as plt
from skimage.registration import phase_cross_correlation
import cv2

def pore_only(Front):
	Front_Pore=np.zeros((Front.shape[0],Front.shape[1]))
	xp=np.arange(0,Front.shape[0])
	for i in range(Front.shape[1]):
		Fp=np.delete(Front[:,i],np.where(np.isnan(Front[:,i])))
		Xp=np.linspace(0,Front.shape[0],len(Fp))
		Front_Pore[:,i]=np.interp(xp,Xp,Fp)
	return Front_Pore

#%% Compute
# import glob, os
# from parse import *
# files_dry=glob.glob(folder_dry+"*."+fmt)
# files_wet=glob.glob(folder_wet+"*."+fmt)
# s=search("{:s}{:d}."+fmt,files_dry[0])

N=np.arange(nmin,nmax,nstep)

fig1,ax1=plt.subplots(int(np.sqrt(len(N)))-1,int(np.sqrt(len(N)))+4,sharey=True,figsize=(10,10)) # rotated interface
fig2,ax2=plt.subplots(1,1) # mean interface
fig3,ax3=plt.subplots(1,1) #  var
fig4,ax4=plt.subplots(1,1) #  hist

MC,VC=[],[]


print("="*(65))

for i,n in enumerate(N):
	print('Frame: ',n)
	
	
	#file1=folder_wet+prefix_wet+'{:04d}.'.format(n)+fmt
	#file2=folder_dry+prefix_dry+'{:04d}.'.format(n)+fmt

#	try:
	#I1=np.float32(cv2.imread(file1,2))/2**16
	#I2=np.float32(cv2.imread(file2,2))/2**16
	
	I1=np.float32(File_wet[args.hdf5name][:,:,n])/2**16
	I2=np.float32(File_dry[args.hdf5name][:,:,n])/2**16
	
	if len(reg)>0:
		# Small range for registration
		I2s=I2[sX[0]:sX[1],sY[0]:sY[1]]
		I1s=I1[sX[0]:sX[1],sY[0]:sY[1]]
		# registration
		shift, error, diffphase = phase_cross_correlation(np.uint8(I1s>th)*255, np.uint8(I2s>th)*255)
		print('Registration - Computed shift (px) :', shift)
		bnd=50
		I1=I1[bnd+int(shift[0]):-bnd+int(shift[0]),bnd+int(shift[1]):-bnd+int(shift[1])]
		I2=I2[bnd:-bnd,bnd:-bnd]

	x,y=np.meshgrid(np.arange(I2.shape[1]),np.arange(I2.shape[0]))
	
	xc=(I2.shape[1])/2
	yc=(I2.shape[0])/2
	
	mask=(x-xc)**2+(y-yc)**2<rad**2
	
	issolid=I2*mask>0.9*np.median(I2[mask])
	isfluid=~issolid*mask
	
	#plt.imshow(issolid)
	
	I=I1*isfluid
	
	#plt.imshow(np.log(I*isfluid))
	
	#plt.figure()
	ax4.hist(I[isfluid].flatten(),1000,density=True,alpha=0.2,color=plt.cm.cool(i/len(N)),label='z={:d}'.format(n))
	
	
	# rescale intensity
	Ir=np.minimum(np.maximum((I-p1)/(p2-p1),0),1)
	
	#plt.figure()
	#ax1[i].imshow(Ir)
	
	#FInd Center of mass
	
	
	if args.rot==0: # automatic rotation calculus
		xm=np.mean(Ir[Ir>0.5]*x[Ir>0.5])/np.mean(Ir[Ir>0.5])-xc
		ym=np.mean(Ir[Ir>0.5]*y[Ir>0.5])/np.mean(Ir[Ir>0.5])-yc
		
		# Rotation to get front straight
		angle=np.arctan(ym/xm)*180/np.pi
	else: # prescribe rotation (in degree)
		angle=np.float16(args.rot)

	If=Ir
	If[~isfluid]=np.nan
	print('Front rotation angle : ', angle)
	from scipy.ndimage import rotate
	
	xi=1+np.arange(Ir.shape[1])
	xi=np.uint32(xi/(Ir.shape[1]/dn))
	x=np.tile(xi.reshape(-1,Ir.shape[1]), (Ir.shape[0],1))
	xr=np.uint32(rotate(x,-angle,reshape=False,cval=np.nan,mode='constant'))
	
	
	# compute binned statistics around the front (0.25-0.75)
	dnfront=np.arange(dn//4,3*dn//4)
	MeanC=np.zeros(len(dnfront))+np.nan
	MedianC=np.zeros(len(dnfront))+np.nan
	varC=np.zeros(len(dnfront))+np.nan
	for j,k in enumerate(dnfront):
		Mask=xr==k
		MeanC[j]=np.nanmean(If[Mask])
		varC[j]=np.nanvar(If[Mask])
		MedianC[j]=np.nanmedian(If[Mask])
	nf=100
	Front=If[If.shape[0]//2-nf:If.shape[0]//2+nf,If.shape[1]//2-nf:If.shape[1]//2+nf]
	
# 		MaskFront=(xr>dn*0.4)&(xr<dn*0.6)
# 		nfront=MaskFront.max()-MaskFront.min()+1
# 		Front=Ir[MaskFront].reshape(nfront,-1)
		
	#Plots
	X=np.arange(len(MeanC))*I.shape[1]/dn
	ax2.plot(X,MeanC,color=plt.cm.cool(i/len(N)),label='z={:d}'.format(n),alpha=0.3)
	ax3.plot(X,varC,color=plt.cm.cool(i/len(N)),label='z={:d}'.format(n))
	ax1.flatten()[i].imshow(pore_only(Front),cmap='inferno')
	ax1.flatten()[i].text(0,0,'{:d}'.format(n))
	ax1.flatten()[i].axis('off')
	
	# Save Z variables
	MC.append(MeanC)
	VC.append(varC)

ax4.set_xlabel('Concentration')
ax4.set_ylabel('Probability')
ax4.set_yscale('log')
ax3.set_ylabel('Variance of C')
ax3.set_xlabel('Distance X accross front')

import os

# Save figures
folder_result=path+'/'+args.SAMPLE+'/'+args.SAMPLE+'_'+args.TYPE+'_'+args.RESOLUTION +'_PostProc/'

try:
	os.mkdir(folder_result)
except:
	print('Result folder exist, erasing...')

fig1.savefig(folder_result+'Front_2D.png')
fig2.savefig(folder_result+'Front_mean.png')
fig3.savefig(folder_result+'Front_variance.png')
fig4.savefig(folder_result+'Front_histogram.png')

print("="*(65))

print('Computing post treatments...')
# Post treatments
# Fit Mean gradient as a function of distance
dx=len(MC[0])//10
#MC=np.array(MC)
x=np.arange(len(MC[0]))*I.shape[1]/dn
grad,amplitude,vc=[],[],[]
#plt.figure()
for i,mc in enumerate(MC):
	left=np.mean(mc[:dx])
	right=np.mean(mc[-dx:])
	amplitude.append(right-left)
	#print(left,right)
	idfit=np.where((mc>right-amplitude[-1]/4)&(mc<left+amplitude[-1]/4))[0]
	if len(idfit)>0:
		p=np.polyfit(x[idfit],mc[idfit],1)
		grad.append(p[0])
		#ax2.plot(x,mc,color=plt.cm.cool(i/len(MC)),alpha=0.2)
		ax2.plot(x[idfit],p[0]*x[idfit]+p[1],'--',color=plt.cm.cool(i/len(MC)))
		vc.append(np.mean(VC[i][idfit]))
	else:
		grad.append(np.nan)
		vc.append(np.nan)
		
scale=np.array(amplitude)/np.array(grad)
ax2.set_xlabel('X (px)')
ax2.set_ylabel('Concentration')

fig,ax=plt.subplots(4,1,sharex=True,figsize=(5,10))
ax[0].plot(N,scale)
ax[0].plot((N+1),10*(N+1)**0.5,'r--')
ax[0].set_ylabel('Mean Front width (px)')


ax[1].plot(N,np.abs(grad))
ax[1].plot((N+1),0.1/(N+1)**0.5,'r--')
ax[1].set_ylabel('Mean Front grad (px)')


ax[2].plot(N,vc)
ax[2].set_ylabel('Max Front Variance C')


ax[3].plot(N,vc/np.abs(grad))
ax[3].set_ylabel('Variance C / gradC')
ax[3].set_xlabel('Distance Z (px)')
ax[3].set_ylim([0,20])

fig.savefig(folder_result+'Front_stats.png')

print('Results saved in folder :'+folder_result)

# save parameters
with open(folder_result+'parameters.txt','w') as data:
	data.write(str(args.__dict__))


print("="*(65))