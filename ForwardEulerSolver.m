%Forward Euler ODE Solver
%
%Brad Yurgens
%MTH 451, Fall 2013
%
%Inputs
%
%t0 - the initial time state of the system, y(t0) = y0
%
%tf - the final time state of the solution
%
%y0 - initial state vector of the solution
%
%f - the differential equation system as a vector valued function
%
%h - the step size of the solution
%
%Outputs
%
%yfinal - the solution at the final time, tf
%
%yn - the set of all state vectors describing the system at each time-step
function [yfinal, yn] = ForwardEulerSolver(t0,tf,y0,f,h)

hold on;

numSteps = (tf - t0)/h;

yn = zeros(size(f,1),numSteps);

for j = 1:size(f,1)

yn(j,1) = y0(j) + h*f{j}(t0,y0);
    
end    
    
for i = 2:numSteps
   
    for j = 1:size(f,1)
    
    yn(j,i) = yn(j,i-1) + h*f{j}(t0 + (i-1)*h,yn(1:size(f,1),i-1));
    
    scatter(t0 + (i-1)*h,yn(j,i),'green');
    
    end
    
end

yfinal = yn(1:size(f,1),numSteps);

hold off;
end