


% Runge-Kutta(RK4) Method
% param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
% - start_t,end_t : start and end limits of time
% - y0 : Initial condition, y(start_t)=y0
% - stepN : stepN * h = (end_t-start_t) , total EulerExp steps performed
% @ return_RK : Matrix of results. return_MT(:,1) is time t,
% return_RK(:,2) is value y at time
function return_RK = Runge_Kutta_Method(f,start_t,end_t,y0,stepN)
    
end