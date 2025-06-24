#17/01/2323 This code comes from command-tag4.py which si in  /LMCG90/TAGSAM/cluster/.  This code was used to simulate the TAGSAN experiment and now it will be modified to run the simulations for te granular bridges, but with non-spherical particles.
#22/02/2023 Now, I will implement a very rustic version of self-gravity.
#27/03/2023 I need to implement a viscous damping in the code so that energy is removed safely.
#04/07/2023 The code needs to bu update so that it can be used to re-start the simulations.
#19/04/2024 There was a problem when the value of the gravitational constant was changed to its actual value after the settling process.  I will make a smal modification so that the force network is set to zero and then recalculated when the change is done.  Also tol has been changed to 1e-4 and dt to 0.0001 s.
#29/11/2024 The only change here is that gravity is not added to the bridge during the pulling phase.  I am doing this to find out to remove the tensile strength added by graivty.  In the original SSMED simulations, cohesive forces were so large that gravity is minimal, but here, that is not the case.  

import os,sys
import os.path
#sys.path.append('/Users/paul/LMGC90/lmgc90_user_2022/build')
#sys.path.append('/Users/eazema/Documents/Programmes/lmgc90_user_2022-Decembre/build_parallele')
sys.path.append('/Users/paul/LMGC90/lmgc90_dev-dev/build')

from pylmgc90.chipy import *

from numpy import *


Initialize()
#############
#############
#############
#############
# Parametres
#############
# ... de la simu 1
nb_steps_1       = 2000000
freq_outFiles_1  = 5000
freq_display_1   = 5000
vis              = 1.0
#
# ... de la simu 2
nb_steps_2       = 20000000
freq_outFiles_2  = 2000
freq_display_2   = 2000
#
# ... du systeme
file = open('parameters.dat','r')
if file.mode=='r':
    param = file.readlines()
    nb_particles=int(param[0])
    Rmin=double(param[1])
    Rmax=double(param[2])
    Rb=double(param[3])
    FC=double(param[4])
    FoV=int(param[5])
file.close()

Rref=0.0225 #This is a reference average radius in m.  I ran simulations with this particle size (4-5 cm) an the time they took was reasonable

# ... du calcul
theta         = 0.5                # Time integrator
freq_detect   = 1                  # Frequence de detection des contacts
dt = 0.001  #Check!!!!  It was 0.0005
#
#       123456789012345678901234567890
type = 'Stored_Delassus_Loops         '
norm = 'QM/16'
tol = 1e-4
relax = 1.0
gs_it1 = 50                    # Nombre d iteration mini
gs_it2 = 4000                   # Nombre d iteration max = gs_it2 * gs_it1
nlgs_3D_SetWithQuickScramble()
nlgs_3D_DiagonalResolution()
#

### Pour les polyedres
PRPRx_LowSizeArrayPolyr(800)
PRPRx_UseCpCundallDetection(800)
RBDY3_NewRotationScheme()
#############
#############
#############
#############
checkDirectories()             # On test si tout les dossiers sont presents
SetDimension(3)                # Set dimension in chipy for dummies  ###
#utilities_EnableLogMes()       # Gestion des messages de log
utilities_DisableLogMes()       # Gestion des messages de log
#############
### model reading ###
utilities_logMes( 'READ BODIES' )
ReadBodies()
utilities_logMes( 'READ BEHAVIOURS' )
ReadBehaviours()
utilities_logMes('LOAD BEHAVIOURS' )
LoadBehaviours()
#utilities_logMes( 'READ MODELS' )
#ReadModels()
#

############################################################
############################################################
### This is to re-start the simulation if needed ###
############################################################
############################################################

lastfn=0
phase=0
lastpf=0
lastpfi=0
lastt=0

path='./lastfile.dat'
check_file = os.path.isfile(path)
print(check_file)

if check_file:
    file = open('lastfile.dat','r')
    if file.mode=='r':
        LF = file.readlines()
        lastfn = int(LF[0])
        phase  = int(LF[1])
        lastpf = float(LF[2])
        lastpfi= float(LF[3])
        lastt  = float(LF[4])
        lastG  = float(LF[5])
    file.close()

############################################################
############################################################

if lastfn>0:
    print ('Simulation re-started')
    #utilities_logMes('LOAD MODELS' )
    #LoadModels()
    utilities_logMes( 'READ INI DOF' )
    ReadIniDof(lastfn)
    utilities_logMes( 'READ DRIVEN DOF' )
    ReadDrivenDof()
    utilities_logMes( 'LOAD TACTORS' )
    LoadTactors()
    utilities_logMes( 'READ INI Vloc Rloc' )
    ReadIniVlocRloc(lastfn)

else:
    print ('New Simulation')
    #utilities_logMes('LOAD MODELS' )
    #LoadModels()
    utilities_logMes( 'READ INI DOF' )
    ReadIniDof()
    utilities_logMes( 'READ DRIVEN DOF' )
    ReadDrivenDof()
    utilities_logMes( 'LOAD TACTORS' )
    LoadTactors()
    utilities_logMes( 'READ INI Vloc Rloc' )
    ReadIniVlocRloc()


#############
### Fin de la lecture du modele ###
#############

#############
### Ecriture du modele dans OUTBOX ###
#############
WriteBodies()
WriteBehaviours()
WriteDrivenDof()
#############
### Fin de l ecriture du modele ###
#############

#############
### On recupere le nombre total de corps
#############
nb_rbdy3 = RBDY3_GetNbRBDY3()
#############
### definition des parametres du calcul ###
#############
utilities_logMes( 'INIT TIME STEPPING' )
TimeEvolution_SetTimeStep(dt)
Integrator_InitTheta(theta)
### post3D ##

#OpenDisplayFiles()
OpenPostproFiles()
setradius=0.01
PT3Dx_SetDisplayRadius(setradius)
### compute masses ###
ComputeMass()



########################
#Global variables to obtain data
########################
pos=[]
vel=[]
forceG=[]
mass=[]

########################
#Global variables for gravity
########################
fpull1=0.0
fpull2=0.0
G=6.6742e-6 #Gravitational constant

########################
#Global variables for viscous drag
########################
Rmean=0.5*(Rmin+Rmax)
fdrag = [[0 for x in range(6)] for y in range(nb_rbdy3+1)]
rad=[Rmean for x in range(nb_rbdy3+1)]
rad[nb_rbdy3]=0.5
rad[nb_rbdy3-1]=0.5

########################
#This is to obtain the indexes for the boulders
########################
Num_Boulder1 = nb_rbdy3-1   #Bottom boulder
Num_Boulder2 = nb_rbdy3     #Top boulder


########################
#This is to open data files
########################
f_strain=open('strain.dat', 'a')


###############################
###############################









########################
#This is to write data from the simulations
########################
def WRDATA(a,b,t,fp):
    i=0
    for i in range(0,nb_rbdy3+1,1):
        pos.append(0.)
        vel.append(0.)
        mass.append(0.)
        forceG.append(0.)
    
    i=0
    for i in range(0,nb_rbdy3+1,1):
        pos[i]=RBDY3_GetBodyVector('Coor_',i)[0:3]
        vel[i]=RBDY3_GetBodyVector('V____',i)[0:3]
        forceG[i]=RBDY3_GetBodyVector('Fext_',i)[0:6]
        mass[i]=RBDY3_GetMass(i)
    
    if a==1:
        f_pos=open('pos-ini.dat', 'w')
        i=0
        for i in range (1,nb_rbdy3+1,1):
            f_pos.write('%e  %e  %e  %e  %e  %e  %e  %e  %e  %e \n'
                        %(pos[i][0],pos[i][1],pos[i][2],vel[i][0],vel[i][1],vel[i][2],
                          forceG[i][0],forceG[i][1],forceG[i][2],mass[i]))
        f_pos.flush()
        f_pos.close()
    
    if b==1:
        Op=(mass[Num_Boulder1]*pos[Num_Boulder1][1]+mass[Num_Boulder2]*pos[Num_Boulder2][1])/(mass[Num_Boulder1]+mass[Num_Boulder2])
        
        i=0
        R_bridge2=0.0
        for i in range(1,nb_rbdy3-1,1):  #It is nb_rbdy3-1 because it only concerns the regolith, not the boulders
            if math.fabs(pos[i][1]-Op)<Rmax:
                dxz2=pos[i][0]*pos[i][0]+pos[i][2]*pos[i][2]
                if (dxz2>R_bridge2):
                    R_bridge2=dxz2
    
        A_bridge=math.pi*R_bridge2
        
        R_bridge=pow(dxz2,0.5)
        f_strain.write('%lf  %lf  %lf  %lf  %e  %e  %e  %lf  %lf\n' %(t,pos[Num_Boulder1][1],pos[Num_Boulder2][1],
                                                              pos[Num_Boulder2][1]-pos[Num_Boulder1][1],
                                                              vel[Num_Boulder1][1],vel[Num_Boulder2][1],
                                                              fp,A_bridge,pow(R_bridge2,0.5)))
    f_strain.flush()





########################
#This is to obtain data from the simulations
########################
def OBDATA():
    i=0
    for i in range(0,nb_rbdy3+1,1):
        pos.append(0.)
        vel.append(0.)
        mass.append(0.)
        forceG.append(0.)

    i=0
    for i in range(0,nb_rbdy3+1,1):
        pos[i]=RBDY3_GetBodyVector('Coor_',i)[0:3]
        vel[i]=RBDY3_GetBodyVector('V____',i)[0:3]
        forceG[i]=RBDY3_GetBodyVector('Fext_',i)[0:6]
        mass[i]=RBDY3_GetMass(i)

########################
#This is to calculate gravity
########################

def LRGR():
    i=0
    for i in range(1,nb_rbdy3+1,1):
        pos[i]=RBDY3_GetBodyVector('Coor_',i)[0:3]
        vel[i]=RBDY3_GetBodyVector('V____',i)[0:3]
    
    d=0
    d2=0
    cosa=0
    cosb=0
    cosc=0
    dx=0
    dy=0
    dz=0
    fgr=0
    fx=0
    fy=0
    fz=0

    
    i=0
    for i in range(0,nb_rbdy3+1,1):
        forceG[i][0]=0.0
        forceG[i][1]=0.0
        forceG[i][2]=0.0
        forceG[i][3]=0.0
        forceG[i][4]=0.0
        forceG[i][5]=0.0
    
    i=0
    j=0
    for i in range (1,nb_rbdy3+1,1):
        for j in range (i+1,nb_rbdy3+1,1):
            dx=pos[j][0]-pos[i][0]
            dy=pos[j][1]-pos[i][1]
            dz=pos[j][2]-pos[i][2]
            d2=dx*dx+dy*dy+dz*dz
            d=pow(d2,0.5)
            cosa=dx/d
            cosb=dy/d
            cosc=dz/d
            fgr=G*mass[i]*mass[j]/d2
            fx=fgr*cosa
            fy=fgr*cosb
            fz=fgr*cosc
            forceG[i][0]+=fx
            forceG[i][1]+=fy
            forceG[i][2]+=fz
            forceG[j][0]+=-fx
            forceG[j][1]+=-fy
            forceG[j][2]+=-fz

    #This chechs if any particle is inside the spheres and moves it out, at the same height.
    i=0
    j=0
    for i in range (1,nb_rbdy3-1,1):
        for j in range (nb_rbdy3-1,nb_rbdy3+1,1):
            dx=pos[j][0]-pos[i][0]
            dy=pos[j][1]-pos[i][1]
            dz=pos[j][2]-pos[i][2]
            d2=dx*dx+dy*dy+dz*dz

            if d2<=Rb*Rb:
                print('Cerca!!!')
                d=pow(d2,0.5)
                cosa=dx/d
                cosb=dy/d
                cosc=dz/d
                print(i,'a x=',pos[i][0],' y=',pos[i][1],' z=',pos[i][2])

                RBDY3_PutBodyVector("Xbeg_",i,[-Rb*cosa,pos[i][1],-Rb*cosc,0.,0.,0.])
                RBDY3_PutBodyVector("Vbeg_",i,[0.,0.,0.,0.,0.,0.])
                print('ca=',cosa,'   cb=',cosb,'  cc=',cosc)
                pos[i]=RBDY3_GetBodyVector('Coor_',i)[0:3]
                print(i,'d x=',pos[i][0],' y=',pos[i][1],' z=',pos[i][2])




########################
#This is the calculate viscous drag
########################
def STDRAG():
    
    i=0
    j=0

    for i in range (1,nb_rbdy3+1,1):
        vptc=RBDY3_GetBodyVector('V____',i)
        fdrag[i][0]=-vis*rad[i]*vptc[0]
        fdrag[i][1]=-vis*rad[i]*vptc[1]
        fdrag[i][2]=-vis*rad[i]*vptc[2]
        fdrag[i][3]=0.0
        fdrag[i][4]=0.0
        fdrag[i][5]=0.0

    fdrag[nb_rbdy3-1][1]=10.0*fdrag[nb_rbdy3-1][1]
    fdrag[nb_rbdy3][1]=10.0*fdrag[nb_rbdy3][1]





########################
#This is the settling routine
########################

def SETTLING():
    
    bulk_behav_SetGravity([0.,0.,0.])
    ts=0.0

    #####################################################
    #Allow DOF=2 in velocity (this is vy) to move freely#
    if FoV==1:
        RBDY3_SetInvisibleVlocyDrivenDof(Num_Boulder1,2)
        RBDY3_SetInvisibleVlocyDrivenDof(Num_Boulder2,2)
    #####################################################
    
    #####################################################
    ### This is needed if re-starting the simulation ###
    if phase==0:
        global lastfn
        if lastfn>0:
            ts=lastt
            print ('ts has been set to ',ts,'s')
    #####################################################

    dKE = 1000.0 #This is any big enough number
    epsi=1e-5
    set_zerov=1
    LRGR() #First call to gravity
    k=0
    kk=0
    tcG=0.0 #Time when G is changed
    tgr1=0.5 #How often gravity is calculated for settling.  Check!!!!
    counter=0
    steps_gr1=int(tgr1/dt)

    #############
    ### DEBUT DU CALCUL 1 : the granular bridge is formed
    #############

    print ('The grains are settling\n')

    for k in range(1, nb_steps_1+1, 1):
        ts += dt
        
        global G
        if G==6.6742e-11:
            LRGR()
            break
        
        KE = postpro_3D_GetKineticEnergy() #Kinetic energy at the beginning
        IncrementStep()
        utilities_logMes( 'COMPUTE Fext' )
        ComputeFext()
        ###############
        
        #Add gravity
        kk+=1
        if kk==steps_gr1:#Gravity is calculated every 100 iterations
            print ('k =',k,'\tts =',round(ts,4),'s \tKE =',format(KE,"0.3e"),'J \tdKE =',format(dKE,"0.3e"),'J')

            print(kk, ' Calculating self-gravity')
            LRGR()
            kk=0
            fpull1=forceG[Num_Boulder1][1]
            fpull2=forceG[Num_Boulder2][1]
    #        print ('fp1=',fpull1,'  fp2=',fpull2)

        STDRAG() #Drag is calculated every time step
        ii=0
        for ii in range(1,nb_rbdy3+1,1): #Gravity and viscous drag added in every time step
            RBDY3_PutBodyVector('Fext_',ii, forceG[ii])
            RBDY3_PutBodyVector('Fext_',ii, fdrag[ii])

        ###############
        utilities_logMes( 'COMPUTE Fint' )
        ComputeBulk()
        utilities_logMes( 'COMPUTE Free Vlocy' )
        ComputeFreeVelocity()
        utilities_logMes( 'SELECT PROX TACTORS' )
        SelectProxTactors()
        RecupRloc()
        ExSolver(type,norm, tol, relax, gs_it1, gs_it2)
        StockRloc()
        utilities_logMes( 'COMPUTE DOF' )
        ComputeDof()
        utilities_logMes( 'UPDATE DOF' )
        UpdateStep()
        ### post3D ###
        WriteOutDof(freq_outFiles_1)
        WriteOutVlocRloc(freq_outFiles_1)
 #       WriteDisplayFiles(freq_display_1)
        WritePostproFiles()
        ### writeout handling ###
        overall_CleanWriteOutFlags()
        
        ### This is to be able to re-start the simulation ###
        counter+=1
        if counter==freq_outFiles_1:
#            global lastfn
            lastfn+=1
            file = open('lastfile.dat','w')
            file.write('%d\n' %lastfn) #Last frame number
            file.write('0\n') #We are still in the settling phase
            file.write('0\n') #Pulling force is zero
            file.write('0\n') #Pulling force increment is zero
            file.write('%f\n' %ts) #Setling time
#            global G
            file.write('%e\n' %G) #Gravitational constant
            counter=0
            file.close()
        
        KE1=postpro_3D_GetKineticEnergy() #Kinetic energy at the beginning
        dKE=math.fabs(KE-KE1)   #Kinectic energy variation
        if dKE<=0.01*epsi and KE<epsi and ts>tcG+1.0:
            print('Settling is finished')
#            global G
            if G>6.6742e-11:
                G=6.6742e-11 #New Gravitational constant
                print('Gravitational Constant has been changed to ',G)
                print('Settling again with new G')
                tcG=ts+5.0
            else:
                LRGR()
                break

                
    WRDATA(1,0,0,0)


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################




##############################################################################

### This is the pulling routine

##############################################################################

def PULLING():

   #####################################
   #####################################

    #This prepares varialbles
    tgr2=5.0 #How often gravity is calculated during the experiment
    steps_gr2=int(tgr2/dt)

    print ('Calculating pulling force')

    fpull1=forceG[Num_Boulder1][1]
    fpull2=forceG[Num_Boulder2][1]
    print('********************')
    print('                    ')
    print ('fp1=',fpull1,'  fp2=',fpull2)
    
    global G
    G=6.6742e-11 #New Gravitational constant

    print('                    ')
    print('                    ')
    print('                    ')
    print('********************')
    LRGR()
    fpull1=forceG[Num_Boulder1][1]
    fpull2=forceG[Num_Boulder2][1]
    print ('fp1=',fpull1,'  fp2=',fpull2,'   G=',G)

    #We know that tensile stress will increase with 1/r^2, then so that simulations do not take too long, I need to increase the pull accordingly.
    #If I consider that the time it took to run a simulation with 4-5cm particles was reasonable, I will take that size as a reference (Rref) and rescale fpinc with respect to that particle size.  That is why I need SzF and why fpinc has SzF as a factor.
    Gfrac=0.02    #Fraction of the gravitational attraction that will be used to pull
    SzF=0.5*(Rmin+Rmax)/Rref    #Particle size factor
    SzF=1.0
    fpinc=(Gfrac/SzF)*0.5*(math.fabs(fpull1)+math.fabs(fpull2))

    print ('fpinc=',fpinc,'N')
    print ('')
    
    #########################################################
    ### If re-started in the settling phase or not at all ###
    #########################################################
    if phase==0:
        print ('Zeroing the system')
        i=0
        for i in range(1,nb_rbdy3+1,1):
            RBDY3_PutBodyVector("Vbeg_",i,[0.,0.,0.,0.,0.,0.])
    #########################################################
    #########################################################
                       
    tc=0.0
    counter=0
    tfp=0.0     #Time to increase fpull
    fpull=0.0   #fpull increment for the entire simulation; this won't change
    tinc=5.0    #Time interval to increase the pulling force
    shift_t=0.0 #This is used in case the simulation is re-started
    shift_kk=0  #This is used in case the simulation is re-started

    ######################################
    ### It re-starting is needed #########
    ######################################
    if phase == 1:
        fpull=lastpf
        tc=lastt
        fpinc=lastpfi
        shift_t=tc%tinc
        if shift_t==0:
            fpull+=fpinc
        shift_kk=int(shift_t/dt)
        tfp=shift_t
        print ('The simulation was re-started and values have')
        print ('been re-established to')
        print ('fpull=',fpull,'N')
        print ('fpinc=',fpinc,'N')
        print ('ts=',tc,'s')
        print ('shift_t=',shift_t)
    
    ####################################
    ###### Pulling speed in y axis #####
    ####################################
    #Make DOF=2 in velocity visible to an imposed velocity
    if FoV==1:
        vpull=2.0e-6
        RBDY3_SetVisibleVlocyDrivenDof(Num_Boulder1, 2)
        RBDY3_SetVisibleVlocyDrivenDof(Num_Boulder2, 2)


    print ('Pulling experiment is starting now...')



    k=0
    kk=shift_kk
    for k in range(1, nb_steps_2+1, 1):

        tc+=dt
        tfp+=dt

        IncrementStep()

        rtc=TimeEvolution_GetTime()

        utilities_logMes( 'COMPUTE Fext' )

        ComputeFext()
        ###############
        #Add gravity
        kk+=1
        if kk==steps_gr2:#Gravity is calculated every tgr2 seconds
            KE = postpro_3D_GetKineticEnergy()
            vb1=RBDY3_GetBodyVector('V____',Num_Boulder1)[0:3]
            print ('k =',k,'\tts =',round(tc,4),'s \tKE =',format(KE,"0.3e"),'J  vb1y=',format(vb1[1],"0.3e"),'m/s')
            
            print(kk, ' Calculating self-gravity')
            #LRGR() #Check!!!
            WRDATA(0,1,tc,fpull)
            kk=0

#        ii=0
#        for ii in range(1,nb_rbdy3+1,1): #Gravity is added in every time step
#            RBDY3_PutBodyVector('Fext_',ii, forceG[ii])

        if FoV==0:
            if tfp>=tinc:    #This increases the pulling force after tinc has passed
                tfp=0.0
                fpull+=fpinc
                print ('fpull has been increased to',format(fpull,"0.3e"),'N')
        
            RBDY3_PutBodyVector('Fext_',Num_Boulder1, [0.0,-fpull,0.0,0.0,0.0,0.0])   #Adding the pulling force in both boulders
            RBDY3_PutBodyVector('Fext_',Num_Boulder2, [0.0,fpull,0.0,0.0,0.0,0.0])
                
        else:
            print('b1=',Num_Boulder1)
            print('b2=',Num_Boulder2)
            RBDY3_SetVlocyDrivenDof(Num_Boulder1,2,-vpull)
            RBDY3_SetVlocyDrivenDof(Num_Boulder2,2,vpull)


        ###############

        utilities_logMes( 'COMPUTE Fint' )
        ComputeBulk()
        utilities_logMes( 'COMPUTE Free Vlocy' )
        ComputeFreeVelocity()
        utilities_logMes( 'SELECT PROX TACTORS' )
        SelectProxTactors()
        RecupRloc()
        nlgs_3D_ExSolver(type,norm, tol, relax, gs_it1, gs_it2)
        StockRloc()
        utilities_logMes( 'COMPUTE DOF' )
        ComputeDof()
        utilities_logMes( 'UPDATE DOF' )
        UpdateStep()

        ### post3D ###
        WriteOutDof(freq_outFiles_2)
        WriteOutVlocRloc(freq_outFiles_2)
#        WriteDisplayFiles(freq_display_2)
        WritePostproFiles()

        ### writeout handling ###
        overall_CleanWriteOutFlags()

        ### This is to be able to re-start the simulation ###
        counter+=1
        if counter==freq_outFiles_2:
            global lastfn
            lastfn+=1
            file = open('lastfile.dat','w')
            file.write('%d\n' %lastfn) #Last frame number
            file.write('1\n') #We are in the pulling phase
            file.write('%e\n' %fpull) #Pulling force
            file.write('%e\n' %fpinc) #Pulling force increment
            file.write('%f\n' %tc) #Pulling time is zero
#            global G
            file.write('%e\n' %G) #Gravitational constant
            counter=0
            file.close()
        
        pyB1=RBDY3_GetBodyVector('Coor_',Num_Boulder1)
        pyB2=RBDY3_GetBodyVector('Coor_',Num_Boulder2)
        if math.fabs(pyB1[1]-pyB2[1])>=1.25:
            break



########################
#THIS IS HOW THE CODE RUNS
########################

if lastfn>0:
    if phase==0:
        OBDATA()
        LRGR()
        STDRAG()
        SETTLING()
        PULLING()
    if phase==1:
        OBDATA()
        LRGR()
        PULLING()

else:
    OBDATA()
    LRGR()
    STDRAG()
    SETTLING()
    PULLING()


#############################
#This is to close the opened files in which I just append data.
f_strain.close()
#############################



WriteOutDof()

WriteOutVlocRloc()

#CloseDisplayFiles()

ClosePostproFiles()

Finalize()


