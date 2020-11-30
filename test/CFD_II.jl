using PyPlot

## Computational Fluid Dynamics 2016 | Coursework 2: Section 3
#NAME: Numerical solution for the Riemann problem
#GOAL: For the linear advection equation; discretised using a the two-step Lax-Friedrichs stencil; appropiate boundary conditions should are applied; then the Riemann problem should be solved to observe the variation of velocity; density; pressure; Mach no. & Entropy/Cv across the computational domain. For the problem is solved for three different mesh sizes.
#Created by: SERGIU PETRE ILIEV. Date: 31.I.2015
#Imperial College London | Aeronautical Engineering | Dr Georgios Papadakis | AE3-414 Computational Fluid Dynamics | Coursework II

# ## Clean & prepare workspace()
# clear               # Clear all variables
# clc                 # Clear the window
# close all           # Close figures
# format long;        # Display more decimal places of the solution

## Initialization of Global Parameters
N=[101,201,301];    # Resolution of the mesh (only 3 mesh density cases will be analysed here)
# N.B. The code has been written so that it can be used for the analysis of a larger number of mesh resolution values - simply include them in the N matrix)
time=0.5;           # Set the time at which the flow parameters should be outputted at the end of the script. The problem will be solved from t=0 to t=time. No need to calculate beyond this point in time.
CFL=1;              # Assign the value of the CFL number [1 is marginally stable, but if the scheme runs correctly it is the most efficient value that can be chosen & still maintain stability]
gamma=1.4;          # Specify the Poisson constant [by default it is set ot 1.4 for air]
L_RHS=2;            # How long the domain is to the right of the burst-disk
L_LHS=2;            # How long the domain is to the left of the burst-disk
L=L_RHS+L_LHS;      # Determine the size of the 1D domain [by default -2 <= x <= 2 implies a domain size of 4]

## Initial Flow conditions - N.B. This could also be written as in the Shock_tube_analytic.m case; where the program asks for pressure & density ratios - the methods are equivalent
# For the driven gas [x .> 0] - at the point immediately to the right of the burst-disk - the script assumes that U_0=U_1
p_r=1;              # Initial gas pressure
rho_r=1;            # Initial gas density N.B. Any gas can be considered, however one should also change the gas constant [gamma] if the program is not applied to air [which is the default]
u_r=0;              # Initial gas velocity

# For the driver gas [x .< 0] - at the point immediately to the left of the burst-disk - the script assumes that U_N+1=U_N
p_l=10;             # Initial gas pressure
rho_l=8;            # Initial gas density
u_l=0;              # Initial gas velocity

## Inialise the matrices that will store the final results for each mesh spacing
p_fin= zeros(size(N,2),maximum(N));     # Initialise pressure array
rho_fin = zeros(size(N,2),maximum(N));  # Initialise density array
u_fin = zeros(size(N,2),maximum(N));    # Initialise velocity array
M_fin= zeros(size(N,2),maximum(N));     # Initialise Mach number array
s_fin= zeros(size(N,2),maximum(N));     # Initialise entropy array
x_fin=zeros(size(N,2),maximum(N));      # Initialise mesh spacing array

## Solving the Riemann problem
# Outer loop to specify the mesh size()
for i=1:1:size(N,2)
    global p, rho, u
    # Initialization of Outer Loop Parameters
    t=0;                # Set the clock back to the beginning of the simulation time
    l=0;                # Set the position tracker back the start of the simulation domain
    NN=2*N[i];          # Since the Lax-Friedrichs scheme is going to be applied, the algorithm requires an efficient way of understanding the half-points i.e. i+1/2 and i-1/2. This can be done switching to a mesh double the size & thinking of each odd point as being a half-point
    dx=L/(NN-1);        # Calculate size of space step depending on number of mesh points
    x=-L_LHS:2*dx:L_RHS;# Generate the 1D mesh in space. Note that it is equispaced with step size dx. Please note that the origin of the system is set at the location of the diaphragm. 

    # Initialize the gas parameters for all points in the mesh; based on the boundary conditions to the RHS & LHS of the diaphragm
    # For the driven gas [x .> 0] - in the domain to the right of the burst-disk
    p[NN/2+1:NN]  = p_r;    # Initial gas pressure
    rho[NN/2+1:NN]= rho_r;  # Initial gas density
    u[NN/2+1:NN]  = u_r;    # Initial gas velocity
    
    # For the driver gas [x .< 0] - in the domain to the left of the burst-disk
    p[1:NN/2]     = p_l;    # Initial gas pressure
    rho[1:NN/2]   = rho_l;  # Initial gas density
    u[1:NN/2]     = u_l;    # Initial gas velocity
       
    c=sqrt(gamma*(p./rho));             # Calculate the speed of sound at each mesh point
    E=(p./rho)*(gamma-1)^(-1)+u.^2/2;   # Calculate the total energy per unit mass of the flow (internal energy + kinetic energy)
    H=E+p./rho;                         # Calculate the total enthalpy per unit mass [total energy + potential energy]
    
    dt=CFL*(dx./max(abs(u)+c));         # Calculate the timestep, as discussed in section 4 of the report

    # Start the Two-Step Lax-Friedrichs solver
    while t<time                        # Begin a while loop iterate the solver for each time step; until the time of interest is reached. N.B. Since t=0 initially; when t<time is no longer true; the last iteration computed was at t=time. 
        l=l+1;                          # Advance the tracker of the position on the 1D mesh()
        
        # Determine the conserved variables for this time instance
        F=[u.*rho p+(u.^2).*rho u.*(p+E.*rho)]';       # Calculate the flux at each point on the mesh()
        U=[rho u.*rho E.*rho]';                        # Calculate flow velocity at each point on the mesh()
        
        # Check if the position tracker is at a half-point | not & apply the appropriate solver to compute U
        if l%1 #if the position tracker is even
            U[1:3,2:NN-2] = (U[1:3,3:NN-1] + U[1:3,1:NN-3])/2  -  dt/(2*dx) * (F[1:3,3:NN-1] - F[1:3,1:NN-3])
        else        #if the position tracker is odd
            U[1:3,3:NN-1] = (U[1:3,4:NN] + U[1:3,2:NN-2])/2  -  dt/(2*dx) * (F[1:3,4:NN] - F[1:3,2:NN-2])
        end
        
        # Find all flow parameters at this timestep across the mesh & store them in their designated arrays
        rho = U[1,1:NN]
        u   = U[2,1:NN]./rho
        E   = U[3,1:NN]./rho
        p   = (E-u.^2/2).*rho*(gamma-1)
        c   = sqrt(gamma*(p./rho))
        
        # Prepare for the next time instance to assure stability
        dt=CFL*(dx./max(abs(u)+c));     # Calculate the next timestep, as discussed in section 4 of the report, note that for stability the timestep will be adaptive as max(abs(u)+c) are likely to change with each iteration
        t=t+dt;                         # Advance the clock one step in time
    end
    
    # Having solved the problem for this specific mesh size; determine & store the flow parameters required for plotting at the time instance of interest i.e. the last time instance where t=time
    rho_fin[i,1:N[i]] = U[1,1:2:NN];                                                                        # Gas density distribution across the mesh at the time of interest
    u_fin[i,1:N[i]]   = U[2,1:2:NN]./rho_fin[i,1:N[i]];                                                     # Gas velocity distribution across the mesh at the time of interest
    p_fin[i,1:N[i]]   = (U[3,1:2:NN]./rho_fin[i,1:N[i]]-u_fin[i,1:N[i]].^2/2).*rho_fin[i,1:N[i]]*(gamma-1); # Gas pressure distribution across the mesh at the time of interest
    M_fin[i,1:N[i]]   = u_fin[i,1:N[i]]./sqrt(gamma*p_fin[i,1:N[i]]./rho_fin[i,1:N[i]]);                    # Local Mach number at each point on the mesh at the time of interest
    s_fin[i,1:N[i]]   = -log(rho_fin[i,1:N[i]])*gamma+log(p_fin[i,1:N[i]]);                                 # Local gas entropy at each point on the mesh at the time of interest
    x_fin[i,1:N[i]]   = x;                                                                                  # Also save the mesh distribution for this iteration
end

## Plot the results for the flow parameters of interest & the time of interest [this section needs to be changed if using a different number of mesh points, by using a simple: for k=1:1:size(N,2) (make, the, figures, using, the, counter, i.e., replacing, the, three, entries, for, each, one, with, one, that, incorporates, the, counter) end, as only the first three will be displayed here - this is because formatting work had to be done in order to differentiate between lines in the report, should it be printed B&W]
#Use the LaTeX interpreter for plotting
set(0,"defaulttextinterpreter','latex") 

figure(1) 
hold all; grid on
title("Velocity Distribution")
plot(x_fin[1,1:N[1]], u_fin[1,1:N[1]],"r-.', 'LineWidth",1.0);      # Plot the numerical solution for N=101
plot(x_fin[2,1:N[2]], u_fin[2,1:N[2]],"b--', 'LineWidth",1.0);      # Plot the numerical solution for N=201
plot(x_fin[3,1:N[3]], u_fin[3,1:N[3]],"k-', 'LineWidth",1.0);       # Plot the numerical solution for N=301
xlabel('x');                                                        # Label the x axis()
ylabel("Velocity");                                                 # Label the y axis()
legend("N=101', 'N=201', 'N=301', 'Location','northeast")           # Specify the legend & its position

figure(2)
hold all; grid on
title("Pressure Distribution")
plot(x_fin[1,1:N[1]], p_fin[1,1:N[1]],"r-.', 'LineWidth",1.0);      # Plot the numerical solution for N=101
plot(x_fin[2,1:N[2]], p_fin[2,1:N[2]],"b--', 'LineWidth",1.0);      # Plot the numerical solution for N=201
plot(x_fin[3,1:N[3]], p_fin[3,1:N[3]],"k-', 'LineWidth",1.0);       # Plot the numerical solution for N=301
xlabel('x');                                                        # Label the x axis()
ylabel("Pressure");                                                 # Label the y axis()
legend("N=101', 'N=201', 'N=301', 'Location','northeast")           # Specify the legend & its position

figure(3)
hold all; grid on
title("Density Distribution")
plot(x_fin[1,1:N[1]], rho_fin[1,1:N[1]],"r-.', 'LineWidth",1.0);    # Plot the numerical solution for N=101
plot(x_fin[2,1:N[2]], rho_fin[2,1:N[2]],"b--', 'LineWidth",1.0);    # Plot the numerical solution for N=201
plot(x_fin[3,1:N[3]], rho_fin[3,1:N[3]],"k-', 'LineWidth",1.0);     # Plot the numerical solution for N=301
xlabel('x');                                                        # Label the x axis()
ylabel("Density");                                                  # Label the y axis()
legend("N=101', 'N=201', 'N=301', 'Location','northeast")           # Specify the legend & its position

figure(4)
hold all; grid on
title("Local Mach number Distribution")
plot(x_fin[1,1:N[1]], M_fin[1,1:N[1]],"r-.', 'LineWidth",1.0);      # Plot the numerical solution for N=101
plot(x_fin[2,1:N[2]], M_fin[2,1:N[2]],"b--', 'LineWidth",1.0);      # Plot the numerical solution for N=201
plot(x_fin[3,1:N[3]], M_fin[3,1:N[3]],"k-', 'LineWidth",1.0);       # Plot the numerical solution for N=301
xlabel('x');                                                        # Label the x axis()
ylabel("Mach number");                                              # Label the y axis()
legend("N=101', 'N=201', 'N=301', 'Location','northeast")           # Specify the legend & its position

figure(5)
hold all; grid on
title("Entropy/Cv Distribution")
plot(x_fin[1,1:N[1]], s_fin[1,1:N[1]],"r-.', 'LineWidth",1.0);      # Plot the numerical solution for N=101
plot(x_fin[2,1:N[2]], s_fin[2,1:N[2]],"b--', 'LineWidth",1.0);      # Plot the numerical solution for N=201
plot(x_fin[3,1:N[3]], s_fin[3,1:N[3]],"k-', 'LineWidth",1.0);       # Plot the numerical solution for N=301
xlabel('x');                                                        # Label the x axis()
ylabel("Entropy/Cv");                                               # Label the y axis()
legend("N=101', 'N=201', 'N=301', 'Location','northeast")           # Specify the legend & its position