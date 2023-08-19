f=@(t,y)y*(2/t+1);

% Euler Explicit Method , using total steps.
% Param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
%
% - start_t,end_t : start and end limits of time
% - start_value : Initial condition, y(start_t)=start_value
% - total_steps : total_steps*h = (end_t-start_t)
% @ return_EES : return Matrix of results. return_EES(:,1) is time,
% return_EES(:,2) is value y at the time.
function return_EES = EulerExpilicit_Steps(f,start_t,end_t,start_value,total_steps)
    h = (end_t-start_t)/total_steps;
    T = linspace(start_t,end_t,total_steps+1);
    Yt = zeros(1,total_steps+1);
    Yt(1) = start_value;
    for j=1:total_steps
        Yt(j+1) = Yt(j)+h*f(T(j),Yt(j));
    end
    return_EES = [T' Yt'];
end
% Euler Explicit Method , using step interval.
% Param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
%
% - start_t,end_t : start and end limits of time
% - start_value : Initial condition, y(start_t)=start_value
% - h : total_steps*h = (end_t-start_t)
% @ return_EES : return Matrix of results. return_EES(:,1) is time,
% return_EES(:,2) is value y at the time.
function return_EEI = EulerExpilicit_Interval(f,start_t,end_t,start_value,h)
    total_steps = ceil((end_t-start_t)/h);
    T = zeros(1,total_steps+1);
    Yt = zeros(1,total_steps+1); 
   % fill time array.
   T(1)=start_t; T(total_steps+1)=end_t;
   for j=1:(total_steps-1)
       T(j+1)=T(j)+h;
   end
   Yt(1) = start_value;
   for j=1:(total_steps-1)
       Yt(j+1) = Yt(j)+h*f(T(j),Yt(j));
   end
   Yt(total_steps+1)=Yt(total_steps)+(T(total_steps+1)-T(total_steps))*f(T(total_steps),Yt(total_steps));
   return_EEI = [T' Yt'];
end

