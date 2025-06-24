# 05-12-2022 This code comes from  gen_sam-cluster7.1.py that was created to simulate the TAGSAM event.
# 06-07-2023 This code has all the routines of gen_sam-bridge2.1.py.  I have made it so that the number of particles, particle size (in cm) and axis ratios are read from a file called particles.ini.  This will facilitate running it in ALPINE as there will be just one master script with different input files.
#23/10/2023 I need to modify the code so that it can run spherical particles too.  The variable PS=0 means spheres, if PS=1, that means polyhedra.

import os,sys
#sys.path.append('/Users/paul/LMGC90/lmgc90_user_2022/build')
sys.path.append('/⁨Users⁩/paulsanchez⁩/⁨Research⁩/⁨3D⁩/⁨lmgc90_dev-dev/build')
import numpy
import math
import random

from pylmgc90 import pre

if not os.path.isdir('./DATBOX'):
    os.mkdir('./DATBOX')



# on se place en 3D
dim = 3

#############
# Parametres
#############
# Particles 0f ~2 cm
#nb_particles     = 5082
#Rmin             = 2.0*0.5*0.01
#Rmax             = 2.0*0.5*0.01
# Particles between 2-3 cm
#nb_particles     = 2373
#Rmin             = 2.0*0.5*0.01
#Rmax             = 3.0*0.5*0.01
# Particles between ~3 cm
#nb_particles     = 2373
#Rmin             = 3.0*0.5*0.01
#Rmax             = 3.0*0.5*0.01
# Particles between 3-4 cm
#nb_particles     = 1287
#Rmin             = 3.0*0.5*0.01
#Rmax             = 4.0*0.5*0.01
# Particles between ~4 cm
#nb_particles     = 1287
#Rmin             = 3.0*0.5*0.01
#Rmax             = 4.0*0.5*0.01
# Particles between 4-5 cm
#nb_particles     = 684
#Rmin             = 4.0*0.5*0.01
#Rmax             = 5.0*0.5*0.01

#Read file for particle properties
file = open('particles.ini','r')
if file.mode=='r':
    particle = file.readlines()
    nb_particles=int(particle[0])   #Number of particles
    Dmin_cm=float(particle[1])      #Minimum  diameter in cm
    Dmax_cm=float(particle[2])      #Maximum diameter in cm
    YR=float(particle[3])           #y/x ratio
    ZR=float(particle[4])           #z/x ratio
    PS=int(particle[5])             #Particle shape, 0=sphere, 1=polyhedron
    FC=float(particle[6])           #Constant cohesive force in N
    FoV=int(particle[7])            #Force[0] or Velocity[1] used to pull boulders apart
file.close()

Rmin             = Dmin_cm*0.5*0.01
Rmax             = Dmax_cm*0.5*0.01

Nb_vertices      = 20
Friction_Part    = 0.6
Cohesion_N       = FC #VdW force for a ~1cm particle
Cohesion_T       = 0.0

#Particles axis relative sizes
#YR               = 1.0
#ZR               = 1.0

#The boulder1 will now be an oblate ellipsoid
Rb=0.5      #1st and 2nd semi-axes of rotation (horizontal plane xy)
Rc=0.5      #3rd semi-axis (vertical axis z)

#Switches
Cs=1        #Cohesion ON (1) or OFF (0)
#############
#############

to_post_process=[]


#Writing parameters to a file so that command.py can read them
file = open('parameters.dat','w')

file.write('%d\n' %nb_particles)
file.write('%lf\n' %Rmin)
file.write('%lf\n' %Rmax)
file.write('%lf\n' %Rb)
file.write('%lf\n' %FC)
file.write('%d\n' %FoV)

file.close()

############## creation des conteneurs
#   * pour les corps
bodies = pre.avatars()
#   * pour les materiaux
mat = pre.materials()
#   * pour les tables de visibilite
svs = pre.see_tables()
#   * pour les lois de contact
tacts = pre.tact_behavs()

############## creations de deux materiaux
#   * un pour les particules
plex = pre.material(name='PLEXx', materialType='RIGID', density=3200)
mat.addMaterial(plex)
#   * un pour les parois
tdur = pre.material(name='TDURx', materialType='RIGID', density=3200)
mat.addMaterial(tdur)

############## on cree un modele de rigide
mod = pre.model(name='rigid', physics='MECAx', element='Rxx3D', dimension=dim)


############## distribution uniforme dans [Rmin, Rmax]
radii=pre.granulo_Uniform(nb_particles, Rmin,Rmax)

############## This is my own routine to deposit the particles
###########All the particles are placed in discs over boulder1
maxlat=Rmax
max_rad=maxlat
r0=Rmax*pow(2.0,0.5)
p=0
q=0
theta_0=0.0
w=1.0

r=[]
theta=[]
rxi=[]
ryi=[]
rzi=[]
coor=numpy.zeros([3*nb_particles])

for p in range(0,nb_particles,1):
    r.append(0)
    theta.append(0)
    rxi.append(0)
    ryi.append(0)
    rzi.append(0)

for p in range(0,nb_particles,1):
    r[p]=r0
    
    
    theta_0=2.0*math.asin(Rmax/r[p])
    theta_0=2.0*math.pi/(numpy.rint(2.0*math.pi/theta_0))
    theta[p]=q*theta_0
    ryi[p]=w*maxlat

    
    if theta[p]>2.001*math.pi-theta_0:
        r0=r[p]+2.2*max_rad
        r[p]=r0
        
        if r[p]+maxlat>Rb:
            r0=Rmax*pow(2.0,0.5)
            w=w+2.0
            if w>5.0:
                print('p=',p)
            ryi[p]=w*maxlat
            r[p]=r[0]
            theta[p]=theta[0]

        
        theta[p]=0.0
        theta_0=0.0
        q=0
    
    q=q+1
    

for p in range(0,nb_particles,1):
    rxi[p]=r[p]*math.cos(theta[p])
    rzi[p]=r[p]*math.sin(theta[p])
#    print (rxi[p],'  ',ryi[p],'  ',rzi[p])

for p in range (0,nb_particles,1):
    coor[3*p]=rxi[p]
    if math.fabs(rxi[p])<0.00001:
        coor[3*p]=0.0
    
    coor[3*p+1]=ryi[p]
    if math.fabs(ryi[p])<0.00001:
        coor[3*p+1]=0.0
    
    coor[3*p+2]=rzi[p]
    if math.fabs(rzi[p])<0.00001:
        coor[3*p+2]=0.0

nb_remaining_particles=nb_particles

############## boucle d'ajout des particules :
if(PS==0):
    print("The particles in this code are spheres.")

    for i in range(0,nb_remaining_particles,1):
        body = pre.rigidSphere(model=mod, material=plex, center=coor[3*i : 3*(i + 1)], color='BLEUx',r=radii[i], number=None)
        # ajout de la sphere dans le conteneur de corps
        bodies += body

if(PS==1):
    print("The particles in this code are polyhedra.")
    
    for i in range(0,nb_remaining_particles,1):
        body = pre.rigidPolyhedron(model=mod, material=plex, center=coor[3*i : 3*(i + 1)],generation_type='regular', color='BLEUx',nb_vertices=Nb_vertices, radius=radii[i],xr=1.,yr=YR,zr=ZR, number=None)
        body.rotate(description='Euler', psi=-random.uniform(0,6.28), center=body.nodes[1].coor)
        body.rotate(description='Euler', theta=-random.uniform(0,6.28), center=body.nodes[1].coor)
        body.rotate(description='Euler', phi=-random.uniform(0,6.28), center=body.nodes[1].coor)
        # ajout de la sphere dans le conteneur de corps
        bodies += body


# On recupere la hauteur maximum
h_max = 0.
for i in range(0,nb_remaining_particles,1):
    h_max_temp = coor[3*i+1] + radii[i]
    if (h_max_temp > h_max):
        h_max = h_max_temp
print ('h_max=',h_max)

########################################
####Boulder creation
########################################


r1=Rb
r2=r1


boulder1=pre.rigidSphere(model=mod, material=plex, center=[0.,0.,0.], color='VERTx',r=Rb, number=None)
if FoV==0:
    boulder1.imposeDrivenDof(component=[1, 3, 4, 5, 6], dofty='vlocy')
if FoV==1:
    boulder1.imposeDrivenDof(component=[1, 2, 3, 4, 5, 6], dofty='vlocy')
else:
    print('Error, FoV is neither 0 nor 1!!!')

import copy
boulder2=copy.deepcopy(boulder1)

bodies+=boulder1
bodies+=boulder2

##############################
#Rotate the boulders so that the smallest tiles do not touch the particles
#boulder1.rotate(description='Euler', psi=0.25*math.pi, axis=[1.,0.0,0.0],center=[0.,0.0,0.0])
#boulder2.rotate(description='Euler', psi=0.25*math.pi, axis=[1.,0.0,0.0],center=[0.,0.0,0.0])
##############################
#On translate tout
#bodies.translate(dx=r2,dy=0.0,dz=r2)
boulder1.translate(dx=0.0,dy=-Rc,dz=0.0)
boulder2.translate(dx=0.0,dy=Rc+h_max,dz=0.0)
bodies.translate(dx=0.0,dy=-0.5*h_max,dz=0.0)
##############################


# gestion des interactions :
#   * declaration des lois


if Cs==0:
    #       - entre particules
    law_part=pre.tact_behav(name='law01',law='RST_CLB',fric=Friction_Part,rstn=0.1,rstt=0.1)
    tacts+=law_part

    #       - avec les boulders
    law_boulder_ptc=pre.tact_behav(name='law02',law='RST_CLB',fric=Friction_Part,rstn=0.1,rstt=0.1)
    tacts+=law_boulder_ptc

    #       - entre les boulders
    law_boulder=pre.tact_behav(name='law03',law='RST_CLB',fric=Friction_Part,rstn=0.1,rstt=0.1)
    tacts+=law_boulder

if Cs==1:
    WTD=0.0005
    #       - entre particules
    law_part=pre.tact_behav(name='law01',law='IQS_WET_DS_CLB',cohn=Cohesion_N,coht=Cohesion_T,Wthk=WTD,dyfr=0.5,stfr=0.6)
    tacts+=law_part

    #       - avec les boulders
    law_boulder_ptc=pre.tact_behav(name='law02',law='IQS_WET_DS_CLB',cohn=1.0*Cohesion_N,coht=Cohesion_T,Wthk=WTD,dyfr=0.5,stfr=0.6)
    tacts+=law_boulder_ptc

    #       - entre les boulders
    law_boulder=pre.tact_behav(name='law03',law='IQS_WET_DS_CLB',cohn=0.0,coht=0.0,Wthk=WTD,dyfr=0.5,stfr=0.6)
    tacts+=law_boulder

##############################
if (PS==0):
    #   * declaration des tables de visibilite
    #       - entre particules
    vis_part = pre.see_table(CorpsCandidat='RBDY3', candidat='SPHER',colorCandidat='BLEUx',
                             behav=law_part, CorpsAntagoniste='RBDY3',
                             antagoniste='SPHER', colorAntagoniste='BLEUx', alert=0.001,halo=1.)
    svs+=vis_part

    #       - entre les particules et les boulders
    vis_boulder_ptc = pre.see_table(CorpsCandidat='RBDY3', candidat='SPHER', colorCandidat='VERTx',
                             behav=law_boulder_ptc, CorpsAntagoniste='RBDY3',
                             antagoniste='SPHER', colorAntagoniste='BLEUx', alert=0.001,halo=1.)
    svs+=vis_boulder_ptc

    #       - entre les boulders
    vis_boulder = pre.see_table(CorpsCandidat='RBDY3', candidat='SPHER', colorCandidat='VERTx',
                                behav=law_boulder, CorpsAntagoniste='RBDY3',
                                antagoniste='SPHER', colorAntagoniste='VERTx', alert=0.001,halo=1.)
    svs+=vis_boulder

if (PS==1):
    #   * declaration des tables de visibilite
    #       - entre particules
    vis_part = pre.see_table(CorpsCandidat='RBDY3', candidat='POLYR',colorCandidat='BLEUx',
                             behav=law_part, CorpsAntagoniste='RBDY3',
                             antagoniste='POLYR', colorAntagoniste='BLEUx', alert=0.001,halo=1.)
    svs+=vis_part
                             
    #       - entre les particules et les boulders
    vis_boulder_ptc = pre.see_table(CorpsCandidat='RBDY3', candidat='SPHER', colorCandidat='VERTx',
                                    behav=law_boulder_ptc, CorpsAntagoniste='RBDY3',
                                    antagoniste='POLYR', colorAntagoniste='BLEUx', alert=0.001,halo=1.)
    svs+=vis_boulder_ptc
     
    #       - entre les boulders
    #vis_boulder = pre.see_table(CorpsCandidat='RBDY3', candidat='SPHER', colorCandidat='VERTx',
    #                            behav=law_boulder, CorpsAntagoniste='RBDY3',
    #                            antagoniste='SPHER', colorAntagoniste='VERTx', alert=0.001,halo=1.)
    #svs+=vis_boulder


# ecriture des fichiers
pre.writeBodies(bodies, chemin='DATBOX/')
pre.writeBulkBehav(mat,dim=3,gravy=[0.,0.,0.],chemin='DATBOX/')
pre.writeTactBehav(tacts, svs, chemin='DATBOX/')
pre.writeDrvDof(bodies, chemin='DATBOX/')
pre.writeDofIni(bodies, chemin='DATBOX/')
pre.writeVlocRlocIni(chemin='DATBOX/')

try:
  pre.visuAvatars(bodies)
except:
  pass

post = pre.postpro_commands()
post.addCommand(pre.postpro_command(name='BODY TRACKING', step=1, rigid_set=to_post_process))
post.addCommand(pre.postpro_command(name='TORQUE EVOLUTION', step=1, rigid_set=to_post_process))
post.addCommand(pre.postpro_command(name='SOLVER INFORMATIONS', step=1))
post.addCommand(pre.postpro_command(name='VIOLATION EVOLUTION', step=1))
pre.writePostpro(commands=post, parts=bodies, path='DATBOX/')
