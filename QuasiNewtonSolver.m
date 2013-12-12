%Newton Solver
%
%
%

function [gn, numIter] = QuasiNewtonSolver(g,J,g0,TOL,MaxIter)

  Jg0 = zeros(size(J,1),size(J,2)); %Initialize container Jacobian and function vector
    
  gk = zeros(size(g,1),1);
    
for numIter = 1:MaxIter
    
    for i = 1:size(J,1)
        
        for j = 1:size(J,1)
            
            Jg0(i,j) = J{i,j}(g0);
        
        end
        
        gk(i) = g{i}(g0);
        
    end
    
    invJg = inv(eye(size(J,1),size(J,2)) - Jg0);
    
    %disp(gk);
    
    %disp(Jg0);
    
    gi = invJg*gk;   %Creat gn+1
    
    if(abs(norm(gi)))< TOL  %Using L-inf norm for 
    
        gn = g0;
        
        return;
        
    else
       
        g0 = gi + g0;
        
    end
end

gn = g0;

end