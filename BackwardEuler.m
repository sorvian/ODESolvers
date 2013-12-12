%Backward Euler

%Inputs
%
%
%
%Outputs
%
%
%

function [yn,ySol] = BackwardEuler(f,J,t0,tf,y0,h)

hold on;

TOL = 10^-4;

maxIter = 10;

numSteps = (tf-t0)/h;

ySol = zeros(size(f,1),numSteps);

ySol(1:size(ySol,1),1) = y0;

tn = t0 + h;

for i = 2:numSteps
    
    g = cell(size(f,1),1);
    
    g0 = zeros(size(g,1),1);
    
    for j = 1:size(f,1) %Generating the vector valued function g for use
                        %with the Quasi Newton Solver
    
        g{j} = @(y) y(j) - h*f{j}(tn,y) - ySol(j,i-1);
        
        g0(j) = ySol(j,i-1) + h*f{j}(tn,ySol(1:size(ySol),i-1));
     
    end
    
    %disp(ySol(1:size(ySol,1),i-1));
    
    yn = QuasiNewtonSolver(g,J,g0,TOL,maxIter);
    
    for j = 1:size(yn,1)
    
        scatter(tn, yn(j), 'red');
    end
        
    tn = t0 + i*h;
    
    ySol(1:size(f,1),i) = yn;
    
    

end

end