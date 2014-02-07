%Runge-Kutta-Fehlberg Fourth and Fifth Order Method for ODEs 
%
%Uses a fifth order method to estimate local truncation error in the fourth
%order method in order to estimate necessary step-size changes, if any based
%upon a devised factor
%
%Inputs
%
%t0 - the initial time state of the interval
%
%tf - the final time state of the interval
%
%y0 - the initial state vector of the system y(t0) = y0
%
%f - the system of differential equation as a vector-valued function
%
%h - the initial time-step
%
%Outputs 
%
%yn - the solution vector at each intervening time-step
%
%numSteps - the number of time-steps taken over the interval


function [yn, numSteps] = RKF45Solver(t0,tf,y0,f,h)

hold  on;

tn = t0;

TOL = 10^-4; %Define a Tolerance, and maximum and minimum time steps for Adaptive Step Sizing

hmin = 0.01;

hmax = 0.25;

numSteps = 0;

while (tn <= tf)
   
   
    fN4 = zeros(1:size(f,1),1); %Generation of RKF4 vector
    
    fN5 = zeros(1:size(f,1),1); %Generation of RKF5 vector
 
    
			%Multi-stage set up of RKF4 and RKF5 methods
    
                                    
           
            for j  = 1:size(f,1)	%k1
			
				y1(j) = h*f{j}(tn,y0);  

			end
    
			for j = 1:size(f,1)		%k2
	
				y2(j) = h*f{j}(tn + 0.25*h,y0 +0.25*y1);
				
			end
			
			for j = 1:size(f,1)		%k3
    
				y3(j) = h*f{j}(tn + 3/8*h,y0 + 3/32*y1 + 9/32*y2);
				
			end
			
			for j = 1:size(f,1)		%k4
    
				y4(j) = h*f{j}(tn + 12/13*h,y0 + 1932/2197*y1 - 7200/2197*y2 + 7296/2197*y3);
				
			end
            
			for j = 1:size(f,1)		%k5
			
				y5(j) = h*f{j}(tn + h, y0 + 439/216*y1 - 8*y2 + 3680/513*y3 - 845/4104*y4);
				
			end
            
            for j = 1:size(f,1)		%k6
			
				y6(j) = h*fn(tn + 0.5*h, y0 - 8/27*y1 + 2*y2 - 3544/2565*y3 + 1859/4104*y4 - 11/40*y5);
			
			end
    
            fN5 = y0 + 16/135*y1 + 6656/12825*y3 + 28561/56430*y4 - 9/50*y5 + 2/55*y6; %Determine y_n+1 for RKF5
            
            fN4 = y0 + 25/216*y1 + 1408/2565*y3 + 2197/4104*y4 - 1/5*y5; %Determine y_n+1 for RKF4
     

 
 eps = abs(max(fN5 - fN4))/h; %Computation of max local truncation error over each component
 
 q = 0.84*(TOL/(eps))^(1/4); %Step-size modification factor
 
 if (q > 4) %Ensure that the convention 0.1 < q < 4 is obeyed
    
     q = 4;
     
 end
 
 if (q < 0.1)
    
     q = 0.1;
     
 end
 
 h = q*h; %adjust step step size
 
 if h > hmax %Ensure that hmin < h < hmax is obeyed
    
     h = hmax;
     
 end
 
 if h < hmin
    
     h = hmin;
     
 end
 
 y0 = fN4; 
    
 tn = tn + h; %Update time value with adjusted step-size
 
  numSteps = numSteps + 1; %increment number of steps;

 for k = 1:size(f,1)
    
    scatter(tn,y0(k),'red');
    
    %scatter(tn,h,'blue');
     
 end
 

 
end

yn = fN4;
end