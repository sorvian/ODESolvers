% PendulumSolver
%
%Function to solve the non-linear, boundary value problem u''(t) = -sin(u(t))
%
%Inputs
%
%Theta_ini - the initial guess for the Jacobian Solver
%
%t0 - the initial time for the solution
%
%T - the final solution time
%
%h - the time step
%
%TOL - the tolerance value for computation 
%
%maxIter - the maximum number of Newton iterations
%
%Outputs
%
%Theta_k - the solution vector for the Jacobian system of the stencil


function Theta_k = PendulumSolver(Theta_ini, t0, T, h, TOL, maxIter)

%Set initial conditions for the Jacobian matrix

Theta_k = Theta_ini;

J = zeros(size(Theta_k,1),size(Theta_k,1));

e = ones(size(Theta_k,1) - 1,1);

%Make every element one row above and below the 

J = J + diag(e,1) + diag(e,-1);


for i = 1:maxIter
   
    %Define the Jacobian via a function handle for the diagonal entries
    %-2 + h^2*cos(Theta_i)
    
    func = @(t) -2 + h.^2*cos(t);
    
    funcTheta_k = func(Theta_k);
    
    J = J + diag(funcTheta_k,0);
    
    %Scale the jacobian to generate the proper stencil
    
    J = 1/h^2*J;
   
    %Solve the system J*delta_k = func(Theta_k)
    
    delta_k = -J/funcTheta_k;
    
    %Check to confirm if the change is within tolerance of L-inf norm
    
    if (abs(norm(delta_k)) <  TOL)
        
        break;
        
    end
    
    Theta_k = Theta_k + delta_k.';
    
end

disp(i);


end