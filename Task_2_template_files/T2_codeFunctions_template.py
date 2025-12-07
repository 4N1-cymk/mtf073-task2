# The header of the corresponding code for Task 1 is relevant also here.
# It is not repeated. Please remind yourself if needed.

# Packages needed
import numpy as np
import matplotlib.pyplot as plt
# Set default font size in plots:
plt.rcParams.update({'font.size': 12})
import os # For saving plots

def createMesh(pointX, pointY,
               mI, mJ, pointXvector, pointYvector):
    ################################
    # DO NOT CHANGE ANYTHING HERE! #
    ################################
    # Only changes arrays in first row of argument list!
    # Sets point coordinates for Task 2 cases.
    for i in range(0, mI):
        for j in range(0, mJ):
            pointX[i,j] = pointXvector[i]
            pointY[i,j] = pointYvector[j]
    
def calcNodePositions(nodeX, nodeY,
                      nI, nJ, pointX, pointY):
    ################################
    # DO NOT CHANGE ANYTHING HERE! #
    ################################
    # Only changes arrays in first row of argument list!
    # Calculates node coordinates.
    # Internal nodes:
    for i in range(0, nI):
        for j in range(0, nJ):
            if i > 0 and i < nI-1:
                nodeX[i,j] = 0.5*(pointX[i,0] + pointX[i-1,0])
            if j > 0 and j < nJ-1:
                nodeY[i,j] = 0.5*(pointY[0,j] + pointY[0,j-1])
    # Boundary nodes:
    nodeX[0,:]  = pointX[0,0]  # Note: corner points needed for contour plot
    nodeY[:,0]  = pointY[0,0]  # Note: corner points needed for contour plot
    nodeX[-1,:] = pointX[-1,0] # Note: corner points needed for contour plot
    nodeY[:,-1] = pointY[0,-1] # Note: corner points needed for contour plot

def calcDistances(dx_PE, dx_WP, dy_PN, dy_SP, dx_we, dy_sn,
                  nI, nJ, nodeX, nodeY, pointX, pointY):
    # Calculate distances in first line of argument list.
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE
    for i in range(1, nI-1):
        for j in range(1, nJ-1):
            dx_PE[i,j] = nodeX[i+1,j] - nodeX[i,j] # ADD CODE HERE x
            dx_WP[i,j] = nodeX[i,j] - nodeX[i-1,j] # ADD CODE HERE x
            dy_PN[i,j] = nodeY[i,j+1] - nodeY[i,j] # ADD CODE HERE x
            dy_SP[i,j] = nodeY[i,j] - nodeY[i,j-1] # ADD CODE HERE x
            dx_we[i,j] = pointX[i,0] - pointX[i-1,0] # ADD CODE HERE x
            dy_sn[i,j] = pointY[0,j] - pointY[0,j-1] # ADD CODE HERE x

def calcInterpolationFactors(fxe, fxw, fyn, fys,
                             nI, nJ, dx_PE, dx_WP, dy_PN, dy_SP, dx_we, dy_sn):
    # Calculate interpolation factors in first row of argument list.
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE
    for i in range(1, nI-1):
        for j in range(1, nJ-1):
            fxe[i,j] = 0.5 * dx_we[i,j] / dx_PE[i,j]  # ADD CODE HERE x
            fxw[i,j] = 0.5 * dx_we[i,j] / dx_WP[i,j] # ADD CODE HERE x
            fyn[i,j] = 0.5 * dy_sn[i,j] / dy_PN[i,j] # ADD CODE HERE x
            fys[i,j] = 0.5 * dy_sn[i,j] / dy_SP[i,j] # ADD CODE HERE x

def initArray(T,
              T_init):
    ################################
    # DO NOT CHANGE ANYTHING HERE! #
    ################################
    # Sets initial default temperature
    # (in all nodes, for contour plot).
    T[:,:] = T_init

def setDirichletBCs(T,
                    nI, nJ, u, v, T_in, T_west, T_east, T_south, T_north):
    # Set Dirichlet boundary conditions
    # (only changes arrays in first row of argument list)
    # Note that a value is needed in all nodes for contour plot
    # Inlets (found by velocity into domain):
    for i in range(nI):
        j = nJ-1
        # ADD CODE HERE x
        # TOP (north) – inlet if v < 0 (downwards into domain)
        if v[i,j] < 0:
            T[i,j] = T_in


        j = 0
        # ADD CODE HERE x SOUTH
        # BOTTOM (south) – inlet if v > 0 (upwards into domain)
        if v[i,j] > 0:
            T[i,j] = T_in
            
        
            
    for j in range(nJ):
        i = nI-1
        # ADD CODE HERE x
        # RIGHT (east) – inlet if u < 0 (into domain)
        if u[i,j] < 0:
            T[i,j] = T_in
        i = 0
        # ADD CODE HERE x
        # LEFT (west) – inlet if u > 0 (into domain)
        if u[i,j] > 0:
            T[i,j] = T_in
            
    # Outlets:
        # Homogeneous Neumann:
        # Set coefficients later, default value already set
    # Walls (found by zero velocity), Dirichlet or initial guess:
    for i in range(nI):
        j = nJ-1
        # ADD CODE HERE x
        # NORTH wall
        if u[i,j] == 0 and v[i,j] == 0:
            T[i,j] = T_north
        j= 0 
        # ADD CODE HERE x
        # SOUTH wall
        if u[i,j] == 0 and v[i,j] == 0:
            T[i,j] = T_south
            
    for j in range(nJ):
        i = nI-1
        # ADD CODE HERE x
        # EAST wall
        if u[i,j] == 0 and v[i,j] == 0:
            T[i,j] = T_east
        i = 0
        # ADD CODE HERE x
        # WEST wall
        if u[i,j] == 0 and v[i,j] == 0:
            T[i,j] = T_west
            
    T[0,0] = T_south
    T[nI-1,nJ-1] = T_north

def calcSourceTerms(Su, Sp,
                    nI, nJ, q_wall, Cp, u, v, dx_we, dy_sn, rho, deltaT, T_o, caseID):
    # Calculate constant source terms
    # (only change arrays in first row of argument list)
    # Keep 'nan' where values are not needed!
    
    # Default values:
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            #no volumetric source
            Su[i,j] = 0 # ADD CODE HERE 
            Sp[i,j] = 0 # ADD CODE HERE
            
    # Heat rate walls (found by zero velocity):
    for i in range(1,nI-1):
        j = nJ-2
        # ADD CODE HERE x
        if q_wall != 0 and abs(u[i, nJ-1]) < 0 and abs(v[i,nJ-1]) < 0:
            Su[i,j] += (q_wall/Cp) * dx_we[i,j]
        
        j = 1
        # ADD CODE HERE x
        if q_wall != 0 and abs(u[i, 0]) < 0 and abs(v[i,0]) < 0:
            Su[i,j] += (q_wall/Cp) * dx_we[i,j]
            
    for j in range(1,nJ-1):
        i = nI-2
        # ADD CODE HERE x
        if q_wall != 0 and abs(u[nI-1, j]) < 0 and abs(v[nI-1,j]) < 0:
            Su[i,j] += (q_wall/Cp) * dy_sn[i,j]
            
        i = 1
        # ADD CODE HERE x
        if q_wall != 0 and abs(u[0, j]) < 0 and abs(v[0,j]) < 0:
            Su[i,j] += (q_wall/Cp) * dy_sn[i,j]
            
    # Time term:
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            #pass # ADD CODE HERE x
            V = dx_we[i,j] * dy_sn[i,j]
            a0P = rho * (V / deltaT)
            Su[i,j] += a0P * T_o[i,j]
            Sp[i,j] -= a0P

def calcD(De, Dw, Dn, Ds,
          gamma, nI, nJ, dx_PE, dx_WP, dy_PN, dy_SP, dx_we, dy_sn):
    # Calculate diffusions conductances in first row of argument list.
    # Note that D is here supposed to include the multiplication with area
    # Only change arrays in first row of argument list!
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE
    for i in range (1,nI-1):
        for j in range(1,nJ-1):
            De[i,j] = gamma * (dy_sn[i,j]/dx_PE[i,j]) # ADD CODE HERE x
            Dw[i,j] = gamma * (dy_sn[i,j]/dx_WP[i,j]) # ADD CODE HERE x
            Dn[i,j] = gamma * (dx_we[i,j]/dy_PN[i,j]) # ADD CODE HERE x
            Ds[i,j] = gamma * (dx_we[i,j]/dy_SP[i,j]) # ADD CODE HERE x

def calcF(Fe, Fw, Fn, Fs,
          rho, nI, nJ, dx_we, dy_sn, fxe, fxw, fyn, fys, u, v):
    # Calculate constant convective (F) coefficients by linear interpolation
    # of velocity in nodes to faces
    # Note that F is here supposed to include the multiplication with area
    # (only changes arrays in first row of argument list)
    # Keep 'nan' where values are not needed!
    # ADD CODE HERE
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            #non-equidistant grid
            u_e = fxe[i,j] * u[i+1,j] + (1 - fxe[i,j]) * u[i,j] #east face normal velocity
            u_w = fxw[i,j] * u[i-1,j] + (1 - fxw[i,j]) * u[i,j] #west face normal velocity
            v_n = fyn[i,j] * v[i,j+1] + (1 - fyn[i,j]) * v[i,j] #north face normal velocity
            v_s = fys[i,j] * v[i,j-1] + (1 - fys[i,j]) * v[i,j] #south face normal velocity
            
            Fe[i,j] = rho * u_e * dy_sn[i,j] # ADD CODE HERE x
            Fw[i,j] = rho * u_w * dy_sn[i,j] # ADD CODE HERE x
            Fn[i,j] = rho * v_n * dx_we[i,j] # ADD CODE HERE x
            Fs[i,j] = rho * v_s * dx_we[i,j] # ADD CODE HERE x

def calcHybridCoeffs(aE, aW, aN, aS, aP,
                     nI, nJ, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs,
                     fxe, fxw, fyn, fys, dy_sn, Sp, u, v,
                     nodeX, nodeY, L, H, caseID):
    # (only changes arrays in first row of argument list)
    # Calculate constant Hybrid scheme coefficients (not taking into account boundary conditions)
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            aE[i,j] = max(0, -Fe[i,j], De[i,j] - (fxe[i,j] * Fe[i,j])) # ADD CODE HERE x
            aW[i,j] = max(0, Fw[i,j], Dw[i,j] + (fxw[i,j] * Fw[i,j])) # ADD CODE HERE x
            aN[i,j] = max(0, -Fn[i,j], Dn[i,j] - (fyn[i,j] * Fn[i,j])) # ADD CODE HERE x
            aS[i,j] = max(0, Fs[i,j], Ds[i,j] + (fys[i,j] * Fs[i,j])) # ADD CODE HERE x
            
    # At outlets (found by velocity out of domain), set homogeneous Neumann
    for j in range(1,nJ-1):
        i = nI-2
        # ADD CODE HERE x
        if u[i+1,j] > 0:
            aE[i,j] = 0
        i = 1
        # ADD CODE HERE x
        if u[i-1,j] < 0:
            aW[i,j] = 0
            
    for i in range(1,nI-1):
        j = nJ-2
        # ADD CODE HERE x
        if v[i,j+1] > 0:
            aN[i,j] = 0
        j = 1
        # ADD CODE HERE x
        if v[i,j-1] < 0:
            aS[i,j] = 0
    
    # (Homogeneous) Neumann walls (found by zero velocity):
    for i in range(1,nI-1):
        j = nJ-2
        # ADD CODE HERE x
        if abs(Fn[i, j]) < 1e-12:
            aN[i,j] = 0
        j = 1
        # ADD CODE HERE
            
    for j in range(1,nJ-1):
        i = nI-2
        # ADD CODE HERE
     
        i = 1
        # ADD CODE HERE x
        if u[i-1,j] == 0:
            aW[i,j] = 0
    
    for i in range(1,nI-1):
        for j in range(1,nJ-1):       
            aP[i,j] = aE[i,j] + aW[i,j] + aN[i,j] + aS[i,j] -Sp[i,j] # ADD CODE HERE

def solveGaussSeidel(phi,
                     nI, nJ, aE, aW, aN, aS, aP, Su, nLinSolIter):
    # Implement the Gauss-Seidel solver for general variable phi,
    # so it can be reused for all variables.
    # Do it in two directions
    # Only change arrays in first row of argument list!
    # ADD CODE HERE
    for linSolIter in range(nLinSolIter):   
        for i in range(1,nI-1):
            for j in range(1,nJ-1):
                #pass # ADD CODE HERE
                #sweep west to east
                phi[i,j] = (aE[i,j]*phi[i+1,j] + aW[i,j]*phi[i-1,j] + aN[i,j]*phi[i,j+1] + aS[i,j]*phi[i,j-1] + Su[i,j]) / aP[i,j]
        for j in range(1,nJ-1):
            for i in range(1,nI-1):
                #pass # ADD CODE HERE
                #sweep south to north
                phi[i,j] = (aE[i,j]*phi[i+1,j] + aW[i,j]*phi[i-1,j] + aN[i,j]*phi[i,j+1] + aS[i,j]*phi[i,j-1] + Su[i,j]) / aP[i,j]
                
def solveTDMA(phi, P, Q,
              nI, nJ, aE, aW, aN, aS, aP, Su, nLinSolIter):
    # Implement the Gauss-Seidel solver for general variable phi,
    # so it can be reused for all variables.
    # Do it in two directions
    # Only change arrays in first row of argument list!
    # ADD CODE HERE
    for linSolIter in range(0,nLinSolIter):
        # March from west to east
        # Sweep from south to north
        for j in range(1,nJ-1):
            # ADD CODE HERE
            i = 1
            a = aP[i,j]
            b = aE[i,j]
            c = aW[i,j]
            d = aN[i,j] * phi[i,j+1] + aS[i,j] * phi[i,j-1] + Su[i,j]
            
            P[i,j] = b / a
            Q[i,j] = (d + c*phi[i-1,j]) / a
            for i in range(2,nI-2):
               # pass # ADD CODE HERE
               a = aP[i,j]
               b = aE[i,j]
               c = aW[i,j]
               d = aN[i,j] * phi[i,j+1] + aS[i,j] * phi[i,j-1] + Su[i,j]
               
               den = a - c*P[i-1,j]
               P[i,j] = b / den
               Q[i,j] = (d + c*Q[i-1,j])/den
               
            #pass# ADD CODE HERE
            i = nI -2
            a = aP[i,j]
            b = aE[i,j]
            c = aW[i,j]
            d = aN[i,j] * phi[i,j+1] + aS[i,j] * phi[i,j-1] + Su[i,j]
            den = a - c*P[i-1,j]
            P[i,j] = 0
            Q[i,j] = (d + c*Q[i-1, j] + b*phi[i+1,j]) / den
            for i in reversed(range(1,nI-1)):
                #pass # ADD CODE HERE
                if i == nI-2:
                    phi[i,j] = Q[i,j]
                else:
                    phi[i,j] = P[i,j]*phi[i+1,j] + Q[i,j]
        # March from north to south
        # Sweep from west to east 
        for i in range(1,nI-1):
            #pass # ADD CODE HERE
            j = 1
            a = aP[i,j]
            b = aN[i,j]
            c = aS[i,j]
            d = aE[i,j] * phi[i+1,j] + aW[i,j] * phi[i-1,j] + Su[i,j]
            
            P[i,j] = b / a
            Q[i,j] = (d + c*phi[i,j-1]) / a
            
            for j in range(2,nJ-2):
                #pass # ADD CODE HERE
                a = aP[i,j]
                b = aN[i,j]
                c = aS[i,j]
                d = aE[i,j] * phi[i+1,j] + aW[i,j] * phi[i-1,j] + Su[i,j]
                
                den = a - c*P[i,j-1]
                P[i,j] = b / den
                Q[i,j] = (d + c*Q[i,j-1]) / den
                
                
            #pass # ADD CODE HERE
            j = nJ - 2
            a = aP[i,j]
            b = aN[i,j]
            c = aS[i,j]
            d = aE[i,j] * phi[i+1,j] + aW[i,j] * phi[i-1,j] + Su[i,j]
            
            den = a - c*P[i,j-1]
            P[i,j] = 0
            Q[i,j] = (d + c*Q[i,j-1] + b*phi[i,j+1]) / den
            
            for j in reversed(range(1,nJ-1)):
                #pass # ADD CODE HERE
                if j == nJ-2:
                    phi[i,j] = Q[i,j]
                    
                else:
                    phi[i,j] = P[i,j]*phi[i,j+1] + Q[i,j]

def correctBoundaries(T,
                      nI, nJ, q_wall, k, dx_PE, dx_WP, dy_PN, dy_SP,
                      u, v, nodeX, nodeY, L, H, caseID):
    # Only change arrays in first row of argument list!
    def is_wall(i,j,u,v,e=1e-12): return (abs(u[i,j]) < e) and (abs(v[i,j]) < e)
    def is_inlet(i,j,u,v,e=1e-12):
        if is_wall(i, j, u, v): return False
        if i == 0: return u[i,j] > e #west:
        if i == nI - 1: return u[i,j] < -e #east:
        if j == 0: return v[i,j] > e #south:
        if j == nJ - 1: return v[i,j] < -e #north:
        return False #interior weird geometry

    def is_outlet(i,j,u,v,e=1e-12):
        if is_wall(i, j, u, v): return False
        if i == 0: return u[i,j] < -e #west:
        if i == nI - 1: return u[i,j] > e #east:
        if j == 0: return v[i,j] < -e #south:
        if j == nJ - 1: return v[i,j] > e #north:
        return False #interior weird geometry
    # Copy T to walls where (non-)homogeneous Neumann is applied
    # Note that specified heat flux is positive INTO computational domain!
    # ADD CODE HERE
    
    # Copy T to outlets (where homogeneous Neumann should always be applied):
    # ADD CODE HERE
    for i in range(1,nI-1):
        j = nJ - 1
        if is_wall(i, j, u, v): T[i,j] = T[i,j-1] # north
        
    for j in range(1,nJ-1):
        i = 0
        if is_wall(i, j, u, v): T[i,j] = T[i+1,j] # west

    
    for i in range(1,nI-1):
    
        j = nJ - 1
        if is_outlet(i, j, u, v): T[i,j] = T[i,j-1] # north
        
    for j in range(1,nJ-1):
        i = 0
        if is_outlet(i, j, u, v): T[i,j] = T[i+1,j] # west
        i = nI - 1
        if is_outlet(i, j, u, v): T[i,j] = T[i-1,j] # east

    # Set cornerpoint values to average of neighbouring boundary points
    T[0,0]   = 0.5*(T[1,0]+T[0,1])     # DO NOT CHANGE
    T[-1,0]  = 0.5*(T[-2,0]+T[-1,1])   # DO NOT CHANGE
    T[0,-1]  = 0.5*(T[1,-1]+T[0,-2])   # DO NOT CHANGE
    T[-1,-1] = 0.5*(T[-2,-1]+T[-1,-2]) # DO NOT CHANGE

def calcNormalizedResiduals(res,
                            nI, nJ, explCorrIter, T,
                            aP, aE, aW, aN, aS, Su, Sp):
    # Compute and print residuals (taking into account normalization):
    # Non-normalized residual:
    r0 = 0.0 # ADD CODE HERE
    for i in range(1,nI-1):
        for j in range(1,nJ-1):
            r = aE[i,j]*T[i+1,j] + aW[i,j]*T[i-1,j] + aN[i,j]*T[i,j+1] + aS[i,j]*T[i,j-1] + Su[i,j]
            l = aP[i,j] * T[i,j]
            r0 += abs(l-r)
    # Append residual at present iteration to list of all residuals, for plotting:
    res.append(r0)
    
    print('iteration: %5d, res = %.5e' % (explCorrIter, res[-1]/res[0]))

def probe(nodeX, nodeY, T, probeX, probeY, method='linear'):
    # method (str): interpolation method ('linear', 'nearest', 'cubic')
    ################################
    # DO NOT CHANGE ANYTHING HERE! #
    ################################
    from scipy.interpolate import griddata
    # Flatten the grid for griddata
    points = np.column_stack((nodeX.ravel(), nodeY.ravel()))
    values = T.ravel()
    # Combine probe coordinates into (M, 2) array
    probes = np.column_stack((probeX, probeY))
    # Perform interpolation
    probe = griddata(points, values, probes, method=method)
    return probe

def createDefaultPlots(
                       nI, nJ, pointX, pointY, nodeX, nodeY,
                       dx_WP, dx_PE, dy_SP, dy_PN, Fe, Fw, Fn, Fs,
                       aE, aW, aN, aS, L, H, T, u, v, k,
                       explCorrIter, res, grid_type, caseID):
    ################################
    # DO NOT CHANGE ANYTHING HERE! #
    ################################
    # (Do not change any input arrays!)
    if not os.path.isdir('Figures'):
        os.makedirs('Figures')

    nan = float("nan")
    
    # Plot mesh
    plt.figure()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Computational mesh \n (Corner nodes only needed for visualization)')
    plt.axis('equal')
    plt.vlines(pointX[:,0],pointY[0,0],pointY[0,-1],colors = 'k',linestyles = 'dashed')
    plt.hlines(pointY[0,:],pointX[0,0],pointX[-1,0],colors = 'k',linestyles = 'dashed')
    plt.plot(nodeX, nodeY, 'ro')
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_mesh.png')
    plt.show()
    
    # Plot velocity vectors
    plt.figure()
    plt.quiver(nodeX.T, nodeY.T, u.T, v.T)
    plt.title('Velocity vectors')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_velocityVectors.png')
    plt.show()
    
    # Plot temperature contour
    plt.figure()
    # plt.contourf(nodeX.T, nodeY.T, T.T)
    tempmap=plt.contourf(nodeX.T,nodeY.T,T.T,cmap='coolwarm',levels=30)
    cbar=plt.colorbar(tempmap)
    cbar.set_label('Temperature $[K]$')
    plt.title('Temperature $[K]$')
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_temperatureDistribution.png')
    plt.show()
    
    # Plot heat flux vectors NORMAL TO WALL boundary face centers ONLY (not in corners)
    # Use temperature gradient just inside domain (note difference to set heat flux)
    qX = np.zeros((nI,nJ))*nan # Array for heat flux in x-direction, in nodes
    qY = np.zeros((nI,nJ))*nan # Array for heat flux in y-direction, in nodes
    for j in range(1,nJ-1):
        i = 0
        if u[i,j] == 0 and v[i,j] == 0:
            dTdx = (T[2,j] - T[1,j]) / dx_PE[1,j]
            qX[i,j] = -k * dTdx # ADD CODE HERE
            qY[i,j] = 0
        i = nI-1
        if u[i,j] == 0 and v[i,j] == 0:
            dTdx = (T[nI-2,j] - T[nI-3,j]) / ( dx_WP[nI-2,j])
            qX[i,j] = -k * dTdx # ADD CODE HERE
            qY[i,j] = 0
    for i in range(1,nI-1):
        j = 0
        if u[i,j] == 0 and v[i,j] == 0:
            dTdy = (T[i,2] - T[i,1]) / ( dy_PN[i,1])
            qX[i,j] = 0
            qY[i,j] = -k * dTdy # ADD CODE HERE
        j = nJ-1
        if u[i,j] == 0 and v[i,j] == 0:
            dTdy = (T[i,nJ-2] - T[i,nJ-3]) / ( dy_SP[i,nJ-2])
            qX[i,j] = 0
            qY[i,j] = -k * dTdy # ADD CODE HERE
            
    plt.figure()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title('Wall-normal heat flux vectors\n (from internal temperature gradient)')
    plt.gca().set_aspect('equal', adjustable='box')
    tempmap=plt.contourf(nodeX.T,nodeY.T,T.T,cmap='coolwarm',levels=30)
    cbar=plt.colorbar(tempmap)
    cbar.set_label('Temperature $[K]$')
    plt.quiver(nodeX, nodeY, qX, qY, color="black")
    plt.xlim(-0.5*L, 3/2*L)
    plt.ylim(-0.5*H, 3/2*H)
    plt.tight_layout()
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_wallHeatFlux.png')
    plt.show()
    
    # Plot residual convergence
    plt.figure()
    plt.title('Residual convergence')
    plt.xlabel('Iterations')
    plt.ylabel('Residual [-]')
    resLength = np.arange(0,len(res),1)
    normalized = [x / res[0] for x in res]
    plt.plot(resLength, normalized)
    plt.grid()
    plt.yscale('log')
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_residualConvergence.png')  
    plt.show()

def createTimeEvolutionPlots(
                             probeX, probeY, probeValues, deltaT, caseID, grid_type):
    # Convert list of arrays to a 2D array: shape (n_steps, n_probes)
    data = np.vstack(probeValues)  # rows = time steps, columns = probe points
    n_steps, n_probes = data.shape
    # Plot evolution for each probe point
    plt.figure()
    for i in range(n_probes):
        plt.plot(range(1, n_steps+1), data[:, i], label=f'Probe {i+1} ({probeX[i]}, {probeY[i]})')
    plt.xlabel('Time Step')
    plt.ylabel('Interpolated Value')
    plt.title('Evolution of Probe Values Over Time')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('Figures/Case_'+str(caseID)+'_'+grid_type+'_timeEvolution.png')
    plt.show()

def createAnimatedPlots(
                       nodeX, nodeY, savedT):
    # Create animated plot
    from matplotlib.animation import FuncAnimation
    fig, ax = plt.subplots()
    # Compute global min and max for consistent color scale
    vmin = min(arr.min() for arr in savedT)
    vmax = max(arr.max() for arr in savedT)
    # Initial contour plot
    tempmap = ax.contourf(nodeX.T, nodeY.T, savedT[0].T,
                          cmap='coolwarm', levels=30, vmin=vmin, vmax=vmax)
    # Add colorbar once
    cbar = fig.colorbar(tempmap, ax=ax)
    cbar.set_label('Temperature [K]')
    # Set static labels and aspect ratio
    ax.set_title('Temperature [K]')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_aspect('equal')
    fig.tight_layout()
    # Update function: redraw contour
    def update(frame):
        ax.clear()  # Clear axis completely
        # Redraw contour for current frame
        ax.contourf(nodeX.T, nodeY.T, savedT[frame].T,
                    cmap='coolwarm', levels=30, vmin=vmin, vmax=vmax)
        # Reset labels and aspect ratio after clearing
        ax.set_title('Temperature [K]')
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_aspect('equal')
        return []  # Nothing to return for blit=False
    # Create animation with looping
    ani = FuncAnimation(fig, update, frames=len(savedT), interval=100,
                        blit=False, repeat=True)
    # plt.show()
    # Save as GIF (works without ffmpeg)
    ani.save('animated_contour.gif', writer='pillow')

def createAdditionalPlots():
    pass

    