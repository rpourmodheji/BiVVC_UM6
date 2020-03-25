import sys
#sys.path.append('/home/likchuan/Research/fenicsheartmodel2')
#sys.path.append('/home/shavik/FenicsHeartCode/artery/lv_fem')

import os as os
from dolfin import *
import numpy as np
from matplotlib import pylab as plt
from petsc4py import PETSc
from forms import Forms
from forms_artery import Forms_artery
from simpleActiveMaterial import SimpleActiveMaterial as Active
from nsolver import NSolver as NSolver
from addfiber_matid import *
import math

parameters["form_compiler"]["quadrature_degree"]=2
parameters["form_compiler"]["representation"] = "uflacs"
#parameters["form_compiler"]["representation"] = "quadrature"
parameters['ghost_mode'] = 'shared_facet'


os.system("rm *.pvd")
os.system("rm *.vtu")
os.system("rm *.pvtu")


# Sequential reading of mesh ###################################################
#mesh = Mesh("../mesh/pmr120_baselinetri.xml")
#facetboundaries = MeshFunction('size_t', mesh, "../mesh/pmr120_baselinetri_facet_region.xml")
#ds = dolfin.ds(subdomain_data = facetboundaries)
#isepiflip=False
#isendoflip=False
#f0, s0, n0 = addfiber_matid(mesh, VectorFunctionSpace(mesh, 'DG', 0), " ", 60, -60, 60, -60, "../", Function(FunctionSpace(mesh,"DG",0)).vector().array()[:], isepiflip, isendoflip)
#File("fiber.pvd") << f0
#File("sheet.pvd") << s0
#File("sheet_normal.pvd") << n0
#topid = 3
#endoid = 2
#epiid = 1
#plot(facetboundaries, interactive=True)
#File("facetboundaries.pvd") << facetboundaries
################################################################################

# Parallel reading of LV mesh ###################################################
meshname = "ellipsoidal"
#meshname = "pmr120"
mesh = Mesh()
f = HDF5File(mpi_comm_world(), "./mesh/"+meshname+".hdf5", 'r')
f.read(mesh, meshname, False)

facetboundaries = MeshFunction("size_t", mesh, 2)
VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=2, quad_scheme="default")
VQuadelem._quad_scheme = 'default'
fiberFS = FunctionSpace(mesh, VQuadelem)

f0 = Function(fiberFS)
s0 = Function(fiberFS)
n0 = Function(fiberFS)

f.read(facetboundaries, meshname+"/"+"facetboundaries")
f.read(f0, meshname+"/"+"eF")
f.read(s0, meshname+"/"+"eS")
f.read(n0, meshname+"/"+"eN")
f.close()

#plot(mesh, interactive=True)
#plot(facetboundaries, interactive = True)

File("facetboundaries.pvd") << facetboundaries
topid = 4
endoid = 1
epiid = 2
########## Reading artery mesh #######################################################
#mesh_art = Mesh("./mesh/artery.xml")
#facetboundaries_art = MeshFunction('size_t', mesh_art, "./mesh/artery_facet_region.xml")
#rightid = 1
#leftid = 2
#innerid = 4
#outerid = 3
##plot(mesh, interactive = True)
##plot(facetboundaries, interactive = True)
##File("mesh.pvd") << mesh
#File("facetboundaries_art.pvd") << facetboundaries_art
meshname_art = "artery"
mesh_art = Mesh()
f_art = HDF5File(mpi_comm_world(), "./mesh/"+meshname_art+".hdf5", 'r')
f_art.read(mesh_art, meshname_art, False)

facetboundaries_art = MeshFunction("size_t", mesh_art, 2)
#matid_art = MeshFunction("size_t", mesh_art)
#f_art.read(matid_art, meshname_art +"/"+"matid_art")
#File("matid_art.pvd") << matid_art

#sub_mesh = SubMesh(mesh_art, matid_art, 0)

VQuadelem_art = VectorElement("Quadrature", mesh_art.ufl_cell(), degree=2, quad_scheme="default")
VQuadelem_art._quad_scheme = 'default'

fiberFS_art = FunctionSpace(mesh_art, VQuadelem_art)

f0_art = Function(fiberFS_art)
s0_art = Function(fiberFS_art)
n0_art = Function(fiberFS_art)

f1_art = Function(fiberFS_art)
s1_art = Function(fiberFS_art)
n1_art = Function(fiberFS_art)

f2_art = Function(fiberFS_art)
s2_art = Function(fiberFS_art)
n2_art = Function(fiberFS_art)

f3_art = Function(fiberFS_art)
s3_art = Function(fiberFS_art)
n3_art = Function(fiberFS_art)

f_art.read(facetboundaries_art, meshname_art +"/"+"facetboundaries")
f_art.read(f0_art, meshname_art+"/"+"eF1")
f_art.read(s0_art, meshname_art+"/"+"eS1")
f_art.read(n0_art, meshname_art+"/"+"eN1")
f_art.read(f1_art, meshname_art+"/"+"eF2")
f_art.read(s1_art, meshname_art+"/"+"eS2")
f_art.read(n1_art, meshname_art+"/"+"eN2")
f_art.read(f2_art, meshname_art+"/"+"eF3")
f_art.read(s2_art, meshname_art+"/"+"eS3")
f_art.read(n2_art, meshname_art+"/"+"eN3")
f_art.read(f3_art, meshname_art+"/"+"eF4")
f_art.read(s3_art, meshname_art+"/"+"eS4")
f_art.read(n3_art, meshname_art+"/"+"eN4")
f_art.close()

File("facetboundaries_art.pvd") << facetboundaries_art
File("mesh_art.pvd") << mesh_art

ds2 = dolfin.ds(subdomain_data = facetboundaries_art)

leftid = 3
rightid = 4
innerid = 1
outerid = 2
##############################################################################
#########################LV model formulation#################################
comm = mesh.mpi_comm()

ls0 = 1.85
N = FacetNormal (mesh)
#Kspring = Constant(100)
X = SpatialCoordinate (mesh)
Press = Expression(("P"), P=0.0, degree=2)
Cavityvol = Expression(("vol"), vol=0.0, degree=2)

Velem = VectorElement("CG", mesh.ufl_cell(), 2, quad_scheme="default")
Velem._quad_scheme = 'default'
Qelem = FiniteElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
Qelem._quad_scheme = 'default'
Relem = FiniteElement("Real", mesh.ufl_cell(), 0, quad_scheme="default")
Relem._quad_scheme = 'default'
Quadelem = FiniteElement("Quadrature", mesh.ufl_cell(), degree=2, quad_scheme="default")
Quadelem._quad_scheme = 'default'

Telem2 = TensorElement("Quadrature", mesh.ufl_cell(), degree=2, shape=2*(3,), quad_scheme='default')
Telem2._quad_scheme = 'default'
for e in Telem2.sub_elements():
	e._quad_scheme = 'default'
Telem4 = TensorElement("Quadrature", mesh.ufl_cell(), degree=2, shape=4*(3,), quad_scheme='default')
Telem4._quad_scheme = 'default'
for e in Telem4.sub_elements():
	e._quad_scheme = 'default'
####### Mixed element for rigid body motion #####################################
VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])
#################################################################################

W = FunctionSpace(mesh, MixedElement([Velem,Qelem,Relem,VRelem]))
Quad = FunctionSpace(mesh, Quadelem)

bctop = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 2), facetboundaries, topid)
bcs = [bctop]

w = Function(W)
dw = TrialFunction(W)
wtest = TestFunction(W)
#du,dp,dpendo, dlsc = TrialFunctions(W)
#(u,p, pendo, lsc) = split(w)
#(v,q, qendo, lsct) = TestFunctions(W)
du,dp,dpendo,dc11 = TrialFunctions(W)
(u,p, pendo,c11) = split(w)
(v,q, qendo,v11) = TestFunctions(W)
#lscprev = Function(Quad)
#lscprev.vector()[:] = ls0
Q = FunctionSpace(mesh,'CG',1)
TF = TensorFunctionSpace(mesh, 'DG', 1)

PK1activestress = Function(Q)
PK1activestress.rename("active stress", "active stress")

t_a = Expression(("t_a"), t_a=0.0, degree=1)
dt = Expression(("dt"), dt=0.0, degree=1)
Tact = Constant(1e5)


params= {"mesh": mesh,
         "facetboundaries": facetboundaries,
         "facet_normal": N,
	 "mixedfunctionspace": W,
	 "mixedfunction": w,
         "displacement_variable": u,
         "pressure_variable": p,
	 "volconst_variable": pendo,
	 "constrained_vol":Cavityvol,
         "endoid": endoid,
	 "fiber": f0,
         "sheet": s0,
         "sheet-normal": n0}

activeparams = {"mesh": mesh,
                "facetboundaries": facetboundaries,
                "facet_normal": N,
                "displacement_variable": u,
                "pressure_variable": p,
                "endoid": endoid,
                "fiber": f0,
                "sheet": s0,
                "sheet-normal": n0,
		"t_a": t_a,
		"dt": dt,
		"Tact": Tact,
		}


uflforms = Forms(params)
activeforms = Active(activeparams)

Fmat = uflforms.Fmat()
Cmat = (Fmat.T*Fmat)
Emat = uflforms.Emat()
J = uflforms.J()

n = J*inv(Fmat.T)*N
dx = dolfin.dx(mesh,metadata = {"integration_order":2})

Ematrix = project(Emat, TF)
Wp = uflforms.PassiveMatSEF()
Wvol = uflforms.V0constrainedE()

Pactive = activeforms.PK1StressTensor()
#lscnew = activeforms.lscnew()


#header_file1 = open("NeoHookeanStress.cpp","r")
#stresscode = header_file1.read()
#StressTensor = dolfin.Expression(stresscode, element= Telem2)
#StressTensor.U = Ematrix
#
#header_file2 = open("NeoHookeanTangent.cpp","r")
#tangentcode = header_file2.read()
#TangentTensor = dolfin.Expression(tangentcode, element= Telem4)
#TangentTensor.U = Ematrix

# Automatic differentiation  #####################################################################################################
F1 = derivative(Wp, w, wtest)
F2 = derivative(Wvol, w, wtest)
#F3 = Kspring*inner(dot(u,n)*n,v)*ds(epiid)
#F3 = Kspring*inner(u,v)*ds(epiid)
F4 = inner(Pactive, grad(v))*dx
L5 = inner(as_vector([c11[0], c11[1], 0.0]), u)*dx + \
	 inner(as_vector([0.0, 0.0, c11[2]]), cross(X, u))*dx + \
	 inner(as_vector([c11[3], 0.0, 0.0]), cross(X, u))*dx + \
	 inner(as_vector([0.0, c11[4], 0.0]), cross(X, u))*dx
F5 = derivative(L5, w, wtest)
Ftotal = F1 + F2 + F4 + F5
#F5 = inner(lsct, lscnew)*dx
#Ftotal = Ftotal + F5

Jac1 = derivative(F1, w, dw)
Jac2 = derivative(F2, w, dw)
#Jac3 = derivative(F3, w, dw)
Jac4 = derivative(F4, w, dw)
Jac5 = derivative(F5, w, dw)
Jac = Jac1 + Jac2 + Jac4 + Jac5
#Jac5 = derivative(F5, w, dw)
#Jac = Jac + Jac5
##################################################################################################################################

# Manual differentiation  ########################################################################################################
#F2 = Press*inner(n, v)*ds(endoid) - Kspring*inner(dot(u,n)*n,v)*ds(epiid) - q*(J-1)*dx - p*J*inner(inv(Fmat.T), grad(v))*dx
#a_metadata = {'quadrature_degree': Telem2.degree(), 'quadrature_scheme': Telem2.quadrature_scheme()}
#F1 = inner(Fmat*StressTensor, grad(v))*dx
#Ftotal = F1 + F2
#
#matrix = 0.5*(grad(du).T*Fmat + Fmat.T*grad(du))
#result = as_tensor(TangentTensor[i,j,k,l]*matrix[k,l],(i,j))
#Jac1 = inner(grad(du)*StressTensor, grad(v))*dx+ grad(v)[i,j]* Fmat[i,k]*result[k,j]*dx
#Jac2 = derivative(F2, w, dw)
#Jac = Jac1 + Jac2
##################################################################################################################################

########################  Artery model Formulation ###################
Cavityart = Expression(("vol_art"), vol_art=0.0, degree=2)
N2 = FacetNormal (mesh_art)
X_art = SpatialCoordinate (mesh_art)

Velem2 = VectorElement("CG", mesh_art.ufl_cell(), 2, quad_scheme="default")
Velem2._quad_scheme = 'default'
Qelem2 = FiniteElement("CG", mesh_art.ufl_cell(), 1, quad_scheme="default")
Qelem2._quad_scheme = 'default'
Relem2 = FiniteElement("Real", mesh_art.ufl_cell(), 0, quad_scheme="default")
Relem2._quad_scheme = 'default'
Quadelem2 = FiniteElement("Quadrature", mesh_art.ufl_cell(), degree=2, quad_scheme="default")
Quadelem2._quad_scheme = 'default'

# Mixed Element for rigid body motion #############################################
VRelem2 = MixedElement([Relem2, Relem2, Relem2, Relem2, Relem2])
###################################################################################
TF_art = TensorFunctionSpace(mesh_art, 'DG', 1)
W2 = FunctionSpace(mesh_art, MixedElement([Velem2, Qelem2, Relem2, VRelem2]))
Quad2 = FunctionSpace(mesh_art, Quadelem2)


#bcleft = DirichletBC(W2.sub(0), Constant(("0.0","0.0","0.0")), facetboundaries_art, leftid)
#bcright = DirichletBC(W2.sub(0), Constant(("0.0","0.0","0.0")), facetboundaries_art, rightid)
bcleft = DirichletBC(W2.sub(0).sub(2), Expression(("0.0"), degree = 2), facetboundaries_art, leftid)
bcright = DirichletBC(W2.sub(0).sub(2), Expression(("0.0"), degree = 2), facetboundaries_art, rightid)
bcs_art = [bcleft, bcright]

#Kspring_art = Constant(("100"))

w2 = Function(W2)
dw2 = TrialFunction(W2)
wtest2 = TestFunction(W2)

du_art, dp_art,dpendo_art, dc = TrialFunctions(W2)
(u_art, p_art, pendo_art,  c1) = split(w2)
(v_art, q_art, qendo_art,  v1) = TestFunctions(W2)

params_art= {"mesh": mesh_art,
         "facetboundaries": facetboundaries_art,
         "facet_normal": N2,
	 "mixedfunctionspace": W2,
	 "mixedfunction": w2,
         "displacement_variable": u_art,
         "pressure_variable": p_art,
	 "volconst_variable": pendo_art,
	 "constrained_vol":Cavityart,
         "innerid": innerid,
	 "fiber1": f0_art,   #45 degree
 	 "fiber2": f1_art,   #135 degree
	 "fiber3": f2_art,   #0 degree, circumferential
	 "fiber4": f3_art,   #90 degree, longitudinal
	}
uflforms2 = Forms_artery(params_art)

Fmat_art = uflforms2.Fmat()
J_art = uflforms2.J()

n2 = J_art*inv(Fmat_art.T)*N2
dx2 = dolfin.dx(mesh_art, metadata = {"integration_order":2})

Wp_art = uflforms2.SEF()*dx2
Wvol_art = uflforms2.V0constrainedE()

# Automatic differentiation  #####################################################################################################
F1_art = derivative(Wp_art, w2, wtest2)
F2_art = derivative(Wvol_art, w2, wtest2)
#F3_art = Kspring_art*inner(u_art,v_art)*ds2(outerid)
# Prevent rigid body motion #################################################################################################
L4_art = inner(as_vector([c1[0], c1[1], 0.0]), u_art)*dx2 + \
	 inner(as_vector([0.0, 0.0, c1[2]]), cross(X_art, u_art))*dx2 + \
	 inner(as_vector([c1[3], 0.0, 0.0]), cross(X_art, u_art))*dx2 + \
	 inner(as_vector([0.0, c1[4], 0.0]), cross(X_art, u_art))*dx2
F4_art = derivative(L4_art, w2, wtest2)
Ftotal_art = F1_art + F2_art + F4_art

Jac1_art = derivative(F1_art, w2, dw2)
Jac2_art = derivative(F2_art, w2, dw2)
#Jac3_art = derivative(F3_art, w2, dw2)
Jac4_art = derivative(F4_art, w2, dw2)
Jac_art = Jac1_art + Jac2_art + Jac4_art
######################################################################################################

solverparams = {"Jacobian": Jac,
                "F": Ftotal,
                "w": w,
                "boundary_conditions": bcs,
		"Type": 1,
		"mesh": mesh
		}
solverparams_art = {"Jacobian": Jac_art,
                "F": Ftotal_art,
                "w": w2,
                "boundary_conditions": bcs_art,
		"Type": 1,
		"mesh": mesh_art
		}

solver= NSolver(solverparams)
solver_art= NSolver(solverparams_art)

########################### Fenics's Newton  #########################################################
Cavityvol.vol = uflforms.cavityvol()
Cavityart.vol_art = uflforms2.cavityvol()
print "unloaded volume LV = ", uflforms.cavityvol()
print "unloaded volume artery = ", uflforms2.cavityvol()*0.001

displacementfile = File("./output/u_disp.pvd")
displacementfile2 = File("./output/u_art.pvd")
activestressfile = File("./output/activestress.pvd")
stressfile = File("./output/art_stress.pvd")

if(MPI.rank(comm) == 0):
	fdataPV = open("PV_.txt", "w", 0)
	fdatals = open("ls.txt","w",0)

# Closed loop cycle
BCL = 800.0
tstep = 0
#Cao = 0.012;
Cven = 0.3;
#Vart0 = 620;#450;
Vven0 = 3200.0;
Rao = 1800.0;
Rven = 2000.0;
Rper = 125000.0;
Rmv = 2000.0;
V_ven = 3700
#V_art = 740;#440
V_LA = 35;

QLVAD = 0.0;
#B_0 = -0.1707;
#B_1 = -0.02177;
#B_2 = 0.0000903;
#omega = 20;
#Rin = 0.10; mmHg*min/L
#Rout = 0.12; mmHg*min/L
Part = 0.0;

cycle = 0.0;
t = 0.0;
tstep = 0.0;
dt.dt = 1.0;
t_lo = 0.0;

LVcav_array = [uflforms.cavityvol()]
Pcav_array = [uflforms.cavitypressure()*0.0075]
Vart_array = [uflforms2.cavityvol()]
Part_array = [uflforms2.cavitypressure()*0.0075]

#Loading Phase
Part = uflforms2.cavitypressure()
V_art = uflforms2.cavityvol()

for p in range(0, 9):
	Cavityart.vol_art += 5000.0
	solver_art.solvenonlinear()
	Part = uflforms2.cavitypressure()
        V_art = uflforms2.cavityvol()
	print "Part = ", Part, " V_art = ", V_art*0.001

# Closed-loop phase
while(cycle < 10):
	#Time varying elastance function for LA and RA
	def et(t, Tmax, tau):
		if (t <= 1.5*Tmax):
			out = 0.5*(math.sin((math.pi/Tmax)*t - math.pi/2) + 1);
      		else:
			out = 0.5*math.exp((-t + (1.5*Tmax))/tau);
		return out

	p_cav = uflforms.cavitypressure()
        V_cav = uflforms.cavityvol()
	Part = uflforms2.cavitypressure()
        V_art = uflforms2.cavityvol()

	tstep = tstep + dt.dt
        cycle = math.floor(tstep/BCL)
	t = tstep - cycle*BCL

	if(t >= 0 and t <4):
		dt.dt = 1.0;
	else:
		dt.dt = 4.0;

	t_a.t_a = t;

#	Part = 1.0/Cao*(V_art - Vart0);
    	Pven = 1.0/Cven*(V_ven - Vven0);
    	PLV = p_cav;

	#### For Calculating P_LA ########################################
        Ees_la = 60;
        A_la = 58.67;
        B_la = 0.049;
        V0_la = 10;
        Tmax_la = 200;
        tau_la = 35;

	if (t < (BCL-200)):
		t_la = t + 200;
        else:
		t_la = t - BCL + 200;

	if(MPI.rank(comm) == 0):
		print "t_LA = ", t_la

        PLA = et(t_la,Tmax_la,tau_la)*Ees_la*(V_LA - V0_la) + (1 - et(t_la,Tmax_la,tau_la))*A_la*(math.exp(B_la*(V_LA - V0_la)) - 1);

	##################################################################################################################################
	if(MPI.rank(comm) == 0):
		print "Cycle number = ", cycle, " cell time = ", t, " tstep = ", tstep, " dt = ", dt.dt
		print >>fdataPV, tstep, p_cav*0.0075 , V_cav, Part*0.0075, V_art*0.001, PLA*0.0075, V_LA

	if(MPI.rank(comm) == 0):
		print "P_ven = ",Pven;
    		print "P_LV = ", PLV;
    		print "P_art = ", Part;
		print "P_LA = ", PLA;

	#### conditions for Valves#######################################
    	if(PLV <= Part):
    	     Qao = 0.0;
    	else:
    	     Qao = 1.0/Rao*(PLV - Part);


    	if(PLV >= PLA):
    	    Qla = 0.0;
    	else:
    	    Qla = 1.0/Rmv*(PLA - PLV);

    	Qper = 1.0/Rper*(Part - Pven);
	Qmv = 1/Rven*(Pven - PLA);
#	QLVAD = 1.0/(Rin + Rout)*(QLVAD*(1 - dt.dt*B_0) + dt.dt*(Part - PLV) - dt.dt*B_2*math.pow(omega,2));      #LVAD flow rate
#	if (cycle >= 8):
#		if (omega == 20.0):
#			QLVAD = (1.0/(Rin + Rout)*(-0.0256*H + 3.2))/60.0;    #20 kRPM
#		if (omega == 22.0):
#			QLVAD = (1.0/(Rin + Rout)*(-0.0236*H + 3.75))/60.0;   #22 kRPM
#		if (omega == 24.0):
#			QLVAD = (1.0/(Rin + Rout)*(-0.0222*H + 4))/60.0;      #24 kRPM
#		if (omega == 26.0):
#			QLVAD = (1.0/(Rin + Rout)*(-0.019*H + 4.4))/60.0;     #26 kRPM
#		if (omega == 28.0):
#			QLVAD = (1.0/(Rin + Rout)*(-0.0176*H + 4.75))/60.0;   #28 kRPM

	if(MPI.rank(comm) == 0):
    		print "Q_LA = ", Qla ;
    		print "Q_ao = ", Qao ;
    		print "Q_per = ", Qper;
		print "Q_mv = ", Qmv ;
		print "QLVAD =", QLVAD;

	V_cav_prev = V_cav
	V_art_prev = V_art
	V_ven_prev = V_ven
	p_cav_prev = p_cav

#	if (cycle >= 8):
#		V_cav = V_cav + dt.dt*(Qla - Qao - QLVAD);   #with LVAD
#		V_art = V_art + dt.dt*(Qao + QLVAD - Qper);  #with LVAD
#	else:
	V_cav = V_cav + dt.dt*(Qla - Qao);
	V_art = V_art + dt.dt*(Qao - Qper)*1000;
    	V_ven = V_ven + dt.dt*(Qper - Qmv);
	V_LA = V_LA + dt.dt*(Qmv - Qla);

	Cavityvol.vol = V_cav
	Cavityart.vol_art = V_art

	if(MPI.rank(comm) == 0):
    		print "V_ven = ", V_ven;
    		print "V_LV = ", V_cav;
    		print "V_art = ", V_art*0.001;
		print "V_LA = ", V_LA;

    	LVcav_array.append(V_cav)
    	Pcav_array.append(p_cav*0.0075)
	Vart_array.append(V_art)
    	Part_array.append(Part*0.0075)

	solver.solvenonlinear()
	solver_art.solvenonlinear()


	ls = sqrt(dot(f0, Cmat*f0))*ls0
	ls1 = project(ls,Q).vector().array()[:]
	eca = project(activeforms.ECa(), Q).vector().array()[:]
	t_r = project(activeforms.tr(), Q).vector().array()[:]
#	ct = project(activeforms.Ct(), Quad).vector().array()[0]
	if(MPI.rank(comm) == 0):
		print >>fdatals, t_a.t_a, min(ls1), max(ls1), min(eca), max(eca), min(t_r), max(t_r)
	#if(MPI.rank(comm) == 0):
	#	print >>fdataPV, t_r
	#lscprev.vector()[:]  = project(activeforms.lsc(), Quad).vector().array()[:]

	#if(MPI.rank(comm) == 0):
	#	print >>output, t_a.t_a, \
	#		        project(ls,Quad).vector().array()[0], \
	#		        project(activeforms.PK1Stress(), Quad).vector().array()[0], \
	#		        project(activeforms.lsc(), Quad).vector().array()[0],  \
	#		        uflforms.cavitypressure()*0.0075, \
	#		        project(activeforms.ftwitch(),Quad).vector().array()[0],\
	#		        project(activeforms.fiso(),Quad).vector().array()[0],\
	#		        project(activeforms.ls_lsc(),Quad).vector().array()[0]



	if(t%4.0 <= 0.0):
        	displacementfile << w.sub(0)
		displacementfile2 << w2.sub(0)
		ls2 = project(ls,Q)
		ls2.rename("sarc", "sarc")
#		sarclenfile << ls2
		PK1activestress.vector()[:] = project(activeforms.PK1Stress(), Q).vector().array()[:]
		activestressfile << PK1activestress
		cauchy_stress = project(uflforms2.PK1stress(),TF_art)
#		print max(cauchy_stress.vector().array()[:])
		cauchy_stress.rename("cauchy_stress", "cauchy_stress")
		stressfile << cauchy_stress



if(MPI.rank(comm) == 0):
	fdataPV.close()
	fdatals.close()
######################################################################################################
