%Rk4 Solver


%Input
%
%
%
%Output
%
%
%


function yn = Rk4Solve(t0,tf,y0,f,h)
    
hold on;    


    tn = t0;
    
    numSteps = (tf - t0)/h;
    
    yn = zeros(size(y0,1),1);
    
    for i = 1:numSteps
      
       for j = 1:size(f,1);                      %Four stage Rk4 set up                              
           
            y1(j) = f{j}(tn,y0);   
            
       end
       
       for j = 1:size(f,1)
                                                   
    
            y2(j) = f{j}(tn + 0.5*h,y0+0.5*h*y1);
            
       end
       
       for j = 1:size(f,1)
    
           y3(j) = f{j}(tn + 0.5*h,y0+0.5*h*y2);
           
       end
    
       for j = 1:size(f,1)
       
            y4(j) = f{j}(tn + h,y0 + h*y3);
            
       end
    
            yn = y0 + h./6.*(y1+2*y2+2*y3+y4); %Determine y_n+1

       for j = 1:size(f,1)
          
           if j == 1
           
           scatter(tn,yn(j),'green');
           
           end
           
           if j == 2
               
               scatter(tn,yn(j),'green');
           end
       end
            
    
    tn = t0 + i*h;
    
    y0 = yn;
    end
    
    hold off;
    
end