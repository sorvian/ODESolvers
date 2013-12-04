








function [yn, iter] = QuasiNewton(f,J,y0,TOL,maxIter)
    
yn = y0;

 for iter =  1:maxIter %Quasi-Newton Rootfinding Method Backward Euler
    
     JN = zeros(size(J,1),size(J,2));
     
     fN = zeros(size(f,1),1);
     
     for j = 1:size(J,1)
       
            for k = 1:size(J,2) %Generate the matrix of Jacobian at y_n
           
               JN(j,k) = J{j,k}(yn); 
                
               %if j == k
                  
                   %JN(j,k) = 1 - JN(j,k); %subtract off the diagonal identity terms
                   
               %end
               
            end
            
             fN(j) = f{j}(yn);
            
     end
     
        vn = -inv(JN)*fN;
       
        if abs(max(vn)) > TOL
           
            yn = yn + vn;
            
        else
           
            break;
            
        end
 end
end