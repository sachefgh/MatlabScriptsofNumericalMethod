f=@(t,y)y*(2/t+1)
E1 = Runge_Kutta(f,1,3,0.3679,3)

% Runge-Kutta Forth-order(RK4)  Method
% param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
% - start_t,end_t : start and end limits of time
% - y0 : Initial condition, y(start_t)=y0
% - stepN : stepN * h = (end_t-start_t) , total EulerExp steps performed
% @ return_RK : Matrix of results. return_MT(:,1) is time t,
% return_RK(:,2) is value y at time
function return_RK = Runge_Kutta(f,start_t,end_t,y0,stepN)
    % Define variables 
    h = (end_t-start_t)/stepN; % step time h
    T = linspace(start_t,end_t,stepN+1); % create&fill T which is step time
    Yt = zeros(1,stepN+1);
    Yt(1) = y0;
    % iteration using RK4 Method
    for i = 1:stepN
        yi05_1 = Yt(i) + 0.5*h*f(T(i),Yt(i));
        yi05_2 = Yt(i) + 0.5*h*f(T(i)+0.5*h,yi05_1);
        yi1_3 = Yt(i) + h*f(T(i)+0.5*h,yi05_2);
        Yt(i+1) = Yt(i) + h*(1/6)*(f(T(i),Yt(i))+2*f(T(i)+0.5*h,yi05_1)+2*f(T(i)+0.5*h,yi05_2)+f(T(i+1),yi1_3));
    end
    return_RK = [T' Yt'];
end