####################################################################################################### 
                                            # MISCELLANEOUS # 
#######################################################################################################   

'''Sampling a Magnetic Field on Grid: Using 2D grids on RZ plane'''
# # generate mesh in R-Z plane at 36 deg
# rrange = np.linspace(1.4, 1.5, 256)
# zrange = np.linspace(-0.2, 0.2, 256)
# R3grid.rzmesh(rrange, zrange, 36.0).savetxt("rzmesh.dat")


# model.load(model = "vessel1_0.05", database = "test")
# # sample magnetic field on mesh
# tasks.magnetic_field("rzmesh.dat")


'''Connection Length: 2D Grid using Integration (NOT APPLICABLE)'''

# get contour of first wall (for systems that have multiple walls or torosurfs) or the innermost vacuum vessel (no need to trace field lines further outside)
# firstwall = boundary.firstwall_rzslice(0) # ANGLE IS IN RADIANS - NEED TO FIND OUT HOW TO CONVERT LATER ON 
# print(type(firstwall)) # Polygon 2D object i.e. contour 


# construct "last closed flux surface" (no need to trace field lines further inside)
# lcfs = PoincareMap.compute(
#     p0 = [????, 0, 0], # p0 is a point estimated from the Poincare plot to be part of the LCFS (USE THE INITIAL PHI FROM POINCARE MAP ALGORITHM) 
#     direction = "fwd", # direction is the direction at which we are tracing field lines toroidally 
#     phiX = 0.0, # phiX is the angle of slice 
#     nsymmetry= 4, 
#     bounded=True, # Whether you want to stop field line tracing at the boundary defined.
#     ).bspline_multifit()
# print(type(lcfs)) # B spline curve object


# generate grid for R-Z slice in divertor region at phi = 0 and exclude unnecessary nodes i.e. outside the boundary or inside the LCFS
# rzmesh = R3grid.rzmesh(
#     rrange = np.linspace(1.27, 1.62, 50), # rrange: Boundary is 5 cm offset in every direction, so we add an extra cm on both ends just to cover it
#     zrange = np.linspace(-0.2, 0.2, 50), # zrange: Again, same thing here 
#     phi = 0.0, #phi section 
#     length_units= 'm',
#     angular_units= 'deg'
#     )

# rzmesh.domain = remove_nodes2d(
#     rzmesh.domain, # Grid 
#     lcfs, # Removal of nodes inside of the grid specified here
#     firstwall # Removal of nodes outside of the grid specified here
#     ) 

# rzmesh.savetxt("../../Data/FLARE_DB/HSX_Test/vessel3_0.15/CL_0.15_rzslice.grid")

# trace field line from grid nodes
# tasks.fieldline_connection(
#     grid = "../../Data/FLARE_DB/HSX_Test/vessel3_0.15/CL_0.15_rzslice.grid", 
#     lcmax = 100, # The maximum connection length: MAIN THING THAT BOTTLENECKS US 
#     output="../../Data/FLARE_DB/HSX_Test/vessel3_0.15/CL_0.15_rzslice.dat", 
#     xfwd = True, # Whether to record the cylindrical coordinates of when the fieldline intersects the boundary in forward direction 
#     xbwd = True, # Same thing but backward direction 
#     alpha = True, # Grazing angle of the fieldline with the boundary 
#     ierr = True # Show errors in calculation 
#     )


'''Connection Length: Field line tracing from LCFS to divertor targets USING INTEGRATION & add output of strike point coordinates on boundary'''
# fieldline_connection(
#     grid = "", #../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSXfluxsurf3d_0.05.grid", 
#     lcmax = 1000, 
#     xfwd = False, 
#     xbwd = False, 
#     ufwd= False, 
#     ubwd= False, 
#     ierr= True,
#     output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/Magnetic_Top_Lc/lc_0.05.dat"
#     ) 
# # THE FWD AND BWD COORDINATES ARE STORED IN THE CURRENT FILE AT WHICH YOU ARE RUNNING THE CODE 

# strike_point_density("lc.dat", dphi=0.4)


'''Field line Diffusion: Field line tracing from LCFS to divertor targets USING INTEGRATION, add output of strike point coordinates on boundary'''

# # Controlling the field line diffusion TO ENABLE IT 
# # control.fieldline.hmax = 0.0175 # limit integration step to ~1 deg
# # control.fieldline.step_type = WHAT IS THE FASTEST? 

# # fieldline_connection(
# #     grid = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/HSXfluxsurf3d_0.05.grid", 
# #     lcmax = 1000, 
# #     xfwd = True, 
# #     xbwd = True, 
# #     ufwd=True, # Needed for strike point density 
# #     ubwd=True, # Needed for strike point density 
# #     ierr= True,
# #     output = "../../Data/FLARE_DB/HSX_Test/vessel1_0.05/lc_0.05.dat"
# #     )