%Forward Euler ODE Solver

function [yfinal, yn] = ForwardEulerSolver(t0,tf,y0,f,h)

numSteps = (tf - t0)/h;

yn = zeros(size(f,1),numSteps);

for j = 1:size(f,1)

yn(j,1) = y0(j) + h*f{j}(t0,y0);
    
end    
    
for i = 2:numSteps
   
    for j = 1:size(f,1)
    
    yn(j,i) = yn(j,i-1) + h*f{j}(t0 + (i-1)*h,yn(1:size(f,1),i-1));
    
    end
    
end

yfinal = yn(1:size(f,1),numSteps);
end