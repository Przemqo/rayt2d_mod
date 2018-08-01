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

#File with layer parameters for unif2
vel_par="vel_par.file"

#Model view *eps
vel_eps="Velocity_Model.eps"
#Output from rayt2d - time model
time_model="nmo_test_model_time.bin"
#Time model view *eps
time_eps="Time_Model.eps"

#Clear previous, if present
rm -f < $vel_model
rm -f < $vel_eps
rm -f < $time_model

#Create dipping model model
#unif2 < $vel_par nx=251 nz=151 dx=20 dz=20 > $vel_model

#Input parameters for Rayt2d
#dt=0.004       #time sample interval in ray tracing             
#nt=401     	#number of time samples in ray tracing      nt=$nt     
                                                                        
fz=0            #first depth sample in velocity                  
nz=151        	# number of depth samples in velocity             
dz=20         	#depth interval in velocity                      
fx=0          	#first lateral sample in velocity                
nx=251        	#number of lateral samples in velocity           
dx=20         	#lateral interval in velocity    

nxs=1			#Number of sources      &
fzs=0         	#Depth of the sources fzs=$fzs
fxs=2500      	#Coordinate of first source fxs=$fxs

fa=-90			#First take-off angle
na=181			#Number of rays
da=1         	#Increment of take-off angle fa=$fa na=$na 


amin=0          #Minimum angle of emergence
amax=180        #Maximum angle of emergence amin=$amin amax=$amax


#smooth2 n1=151 n2=251 r1=5 r2=5 < $vel_model > temp.bin
#mv temp.bin $vel_model

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
rayt2d_mod fz=$fz nz=$nz dz=$dz fx=$fx nx=$nx dx=$dx \
nxs=$nxs fa=$fa na=$na da=$da fxs=$fxs fzs=$fzs amin=$amin amax=$amax \
vfile=$vel_model tfile=$time_model refl="0,720;5000,720" #refl="0,1300;5000,1100" 

#Calculated traveltime cube PS
pscube < $time_model n1=$nz d1=$dz f1=$fz label1="$labelz" \
    n2=$nx d2=$dx f2=$fx label2="$labelx" \
    n3=$nxs d3=$d3 label3="$labely" \
	ybox=3 hbox=4 \
	bclip=2 d1num=500 d2num=1000 d3num=5 \
	title="Traveltime Tables"  > $time_eps


#b2a < $time_model  n1=1 format=1 > nmo_test_model_time.txt
#head -$nz nmo_test_model_time.txt > nmo_test_model_0offset_traveltimes.txt

#Extract traveltimes for the first shot (nx x nz)
#head -37901 nmo_test_model_time.txt > nmo_test_model_traveltimes.txt
#a2b < nmo_test_model_traveltimes.txt n1=1 > nmo_times_shot_1.bin

#View traveltimes from 1st shot
#ximage < nmo_times_shot_1.bin n1=151 n2=251 title="OBS 3007" width=768 height=244 &


