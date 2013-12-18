#Kannst Du die Sim noch mal mit d=0.5...1mm und mu-p 0.9 und mur 0.5 und dt 1e-6
variable scale   equal	 2
variable g	 equal 	 9.81
variable dt	 equal	 1e-6*${scale}

variable factor	 equal	 1/${dt}
variable t1 equal round(0.01*${factor}) #1.2
variable t11 equal round(0.2*${factor}) #0.4

variable t2 equal round(0.9*${factor}) #0.9 .. 1sek Ausflusszeit

#ff 03..06
#frf 02 05

variable 	ff equal   0.15  			# coeff. static friction between particles 
variable 	frf equal  0.2   			# coeff. rolling friction between particles 

variable	wallfric equal 0.1 			# wall friction 	
 
variable 	cor equal 0.735				# coefficient of restitution
variable	N equal    200			# Anzahl Partikel 160000
variable 	dens equal 1040				# Dichte im Labor ermittelt, passt zu HE
variable	poiss equal 0.22
variable	shearmod equal 1e8
variable	youngmod equal 1e7 #${shearmod}*(2.0+2.0*${poiss})
#variable	CED equal 1.3e3 #kpa=J/m³ 1.3kpa by Mitchel et all referenced by He page 68

#remove scaling factor
variable r1 equal 0.00085
variable r2 equal 0.0009
variable r3 equal 0.00095

atom_style	granular
atom_modify	map array
boundary	f f f
newton		off
echo		both

processors 1 1 8    #changed by srl
###################################################################################################
# Simulationsgebiet fuer DEM Partikel
###################################################################################################
communicate	single vel yes
units		si
region		reg block 0 0.3 0 0.02 0 1.0 units box
create_box	2 reg

neighbor	0.0018 bin
neigh_modify	delay 0

#Material properties required for new pair styles

fix 		m1 all property/global youngsModulus peratomtype ${youngmod} ${youngmod}
fix 		m2 all property/global poissonsRatio peratomtype ${poiss} ${poiss}
fix 		m3 all property/global coefficientRestitution peratomtypepair 2 ${cor} ${cor} ${cor} ${cor}
fix 		m4 all property/global coefficientFriction peratomtypepair 2 ${ff} ${wallfric} ${wallfric} ${wallfric}
fix 		m5 all property/global coefficientRollingFriction peratomtypepair 2 ${frf} ${frf} ${frf} ${frf}
fix 		m6 all property/global coefficientRollingViscousDamping peratomtypepair 2 0 0 0 0
#fix 		m7 all property/global cohesionEnergyDensity peratomtypepair 2 ${CED} 0 0 0 #nur cohesion innerhlab des SG ?

#New pair style
pair_style 	gran/hertz/history rolling_friction epsd2 cohesion off #sjkr
pair_coeff	* *

timestep	${dt}

fix		gravi all gravity $g vector 0.0 0.0 -1.0

#granular walls
#fix 		inner all mesh/surface file inner_trans.stl type 2 move -8.782e-3 0 -10e-3
#fix		wall all wall/gran/hertz/history mesh n_meshes 1 meshes inner rolling_friction epsd2 cohesion off #sjkr

####################################################################################################
# DEM Begrenzung durch Waende, primitive type bedeutet ohne besondere Form
####################################################################################################
fix		wall_x1 all wall/gran/hertz/history primitive type 2 xplane 0
fix		wall_x2 all wall/gran/hertz/history primitive type 2 xplane 0.3

fix		wall_y1 all wall/gran/hertz/history primitive type 2 yplane 0
fix		wall_y2 all wall/gran/hertz/history primitive type 2 yplane 0.02

fix		wall_z1 all wall/gran/hertz/history primitive type 2 zplane 0

#fix		lid all wall/gran/hertz/history primitive type 2 zplane 0.079 #(150/2+8/2)


#region and insertion
fix		pts1 all particletemplate/sphere 1 atom_type 1 density constant ${dens} radius constant ${r1}
fix		pts2 all particletemplate/sphere 1 atom_type 1 density constant ${dens} radius constant ${r2}
fix		pts3 all particletemplate/sphere 1 atom_type 1 density constant ${dens} radius constant ${r3}


fix		pdd1 all particledistribution/discrete 1. 3 pts1 1./3. pts2 1./3. pts3 1./3.

#####################################################################################################
# Partikelgenerierung in nachfolgend definiertem Gebiet mit entsprechender Rate
#####################################################################################################
#group		nve_group region reg
region		gen block 0 0.3 0 0.02 0.6 0.8 units box 

#particle insertion
fix ins all insert/rate/region seed 10017 distributiontemplate pdd1 nparticles $N particlerate 1000000 overlapcheck yes vel constant 0 0 -2 insert_every 5000 verbose no region gen 

#apply nve integration to all particles that are inserted as single particles
fix		integr all nve/sphere

#output settings, include total thermal energy
compute		1 all erotate/sphere
fix 		ts_check all check/timestep/gran 1000 0.1 0.1
thermo_style	custom step atoms cpu tpcpu
thermo		1000
thermo_modify	lost ignore norm no
compute_modify	thermo_temp dynamic yes


#insert the first particles so that dump is not empty
shell mkdir post_$g_${ff}_${frf}

fix color all property/atom Color scalar yes no no 0
run 0	#apply

restart 10000 Restart/restart*.bin
dump		dmp1 all custom   1000 post_$g_${ff}_${frf}/dump_*.liggghts id x y z vx vy vz radius f_color
run ${t1}  #every 10000 "loadbalance nlocal/simple ntry 10"
run ${t11} #every 10000 "loadbalance nlocal/simple ntry 10"
unfix ins
#create new particle-property color


#variable c equal 0
#label colorloop
#variable i loop 0 8
# variable zlow  equal 0.079+$i*0.008875
# variable zhigh equal 0.079+$i*0.008875+0.008875
# variable color equal $i
# region R$i block  -0.001 0.101 -1e-16 0.0151 ${zlow} ${zhigh} units box
# set region R$i property/atom Color ${color}
#next	    i
#jump	    in.lmp colorloop

#dump		dmp all custom   1000 post_$g_${ff}_${frf}/dump_*.bin id x y z vx vy vz radius f_color
#write_restart phase1.restart

#unfix lid
dump		dmp all custom   1000 post_$g_${ff}_${frf}/dump_*.liggghts id x y z vx vy vz radius f_color
run 1
#cfd coupling
fix		cfd all couple/cfd couple_every 20 mpi # alle 100 DEM-schritte (muss in couplingProperties CFD geändert werden)
fix		cfd2 all couple/cfd/force		#kräfte aus dem CFD holen

#variable       vx equal vx[1]
#variable       vy equal vy[1]
#variable       vz equal vz[1]
#variable       time equal step*dt
#fix            extra all print 10 "${time} ${vx} ${vy} ${vz}" file ../DEM/post/velocity.txt title "%" screen no



#insert the first particles so that dump is not empty
run		2

#restart 10000 step*.restart
#run ${t2} #every 100000 "loadbalance nlocal/simple ntry 10"

