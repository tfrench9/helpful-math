function [ t, y ]  = rk4( functionString, step, tFinal, initialValue )
% RK4 - Applies Simple Runge-Kutta 4 method for solving ODE.
%   Y = rk4(FUNCTION(t,y), STEP_SIZE, T_FINALS, INITIAL VALUE) Uses step 
%   size of the first input and computes from 0 to T_final.
%
%       dy/dt = F(t,y)
%
%   ex: y = rk4('sin(t*y) + exp(2*t)', 0.01, 5, 0)
%
%       dy/dt = sin(t*y) + exp(2*t)
%
%   this would evaluate this function from t = 1 to 5 using the
%   rk4 method with step size of 0.01
%--------------------------------------------------------------------------
% 
    Function = @(t,y) (eval(functionString));
    t = 0:step:tFinal;
    y = zeros(1, length(t));
    y(1) = initialValue;
   
    for i = 2:numel(t)
        k1 = step * Function( t(i-1) , y(i-1) );
        k2 = step * Function( t(i-1) + step/2, y(i-1) + k1/2 );
        k3 = step * Function( t(i-1) + step/2, y(i-1) + k2/2 );
        k4 = step * Function( t(i-1) + step, y(i-1) + k3);
        
        y(i) = y(i-1) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

end
