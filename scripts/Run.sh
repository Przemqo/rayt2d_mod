#!/bin/bash
#################################################
##############FOR OBS DATA#######################
#################################################
#Converting .tomo velocity format to 2d grid
vtomo2xyz < gr_nmo_test4.tomo > gr_nmo_test4_tomo.txt 

#Extracting velocities from all gridpoints 
awk -F "\"* \"*" '{print $3*1000}' gr_nmo_test4_tomo.txt > gr_nmo_test4_tomo_vmodel.txt 

#Creating velocity model
a2b < gr_nmo_test4_tomo_vmodel.txt n1=1 > nmo_test_model.bin 
#################################################


#Input vel model name
vel_model="nmo_test_model.bin"
#Model view *eps
vel_eps="Velocity_Model.eps"
#Output from rayt2d - time model
time_model="nmo_test_model_time.bin"
#Time model view *eps
time_eps="Time_Model.eps"
#Traveltimes from reflections to surface
traveltimes="traveltimes.txt"

#Clear previous, if present
rm -f < $vel_model
rm -f < $vel_eps
rm -f < $time_model
rm -f < $traveltimes


#Input parameters for Rayt2d
#dt=0.004       #time sample interval in ray tracing             
#nt=401     	#number of time samples in ray tracing     
                                                                        
fz=0            #first depth sample in velocity                  
nz=151        	# number of depth samples in velocity             
dz=20         	#depth interval in velocity                      
fx=0          	#first lateral sample in velocity                
nx=251        	#number of lateral samples in velocity           
dx=20         	#lateral interval in velocity    

nxs=1			#Number of sources      
fzs=0           #Depth of the sources 
fxs=2500      	#Coordinate of first source 

fa=-90			#First take-off angle
na=181			#Number of rays
da=1        	#Increment of take-off angle fa=$fa na=$na 


amin=0          #Minimum angle of emergence
amax=180        #Maximum angle of emergence amin=$amin amax=$amax


#Labels for model ps files
labelz="Depth (m)"
labelx="Distance (m)"
labely="Shot number"

#Velocity model PS file
psimage < $vel_model  style=seismic \
	n1=$nz d1=$dz f1=$fz grid1=dot label1="$labelz" \
	n2=$nx d2=$dx f2=$fx grid2=dot label2="$labelx" \
	wbox=6 hbox=4 ybox=4 \
	d1num=1000 d2num=1000 \
	title="Velocity"  > $vel_eps



# use rayt2d to generate traveltime tables from model
# In refl, You need to specify all the reflectors for which to
# calculate ray paths to and from (they have to be present in velocity model as well)
# They can be just 2 x,z pairs or arbitrary number, dipping or not
# Each reflector as a separate "refl" parameter

rayt2d_mod fz=$fz nz=$nz dz=$dz fx=$fx nx=$nx dx=$dx \
	nxs=$nxs fa=$fa na=$na da=$da fxs=$fxs fzs=$fzs amin=$amin amax=$amax \
	vfile=$vel_model tfile=$time_model rtfile=$traveltimes \
	refl="0,720;5000,720" #refl="0,1300;5000,1100" 

#Calculated traveltime cube PS
pscube < $time_model n1=$nz d1=$dz f1=$fz label1="$labelz" \
    n2=$nx d2=$dx f2=$fx label2="$labelx" \
    n3=$nxs d3=$d3 label3="$labely" \
	ybox=3 hbox=4 \
	bclip=2 d1num=500 d2num=1000 d3num=5 \
	title="Traveltime Tables"  > $time_eps



