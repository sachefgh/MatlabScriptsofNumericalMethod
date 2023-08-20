f=@(t,y)(exp(2*t)*(4*y+6*t-3))/( ((2*y+3*t)^2)+2*exp(2*t) ) %;    % This is function handle of f(t,y)
E1 = ModifiedTraezodial(f,0,2,sqrt(0.5),400) %;
% plot
plot(E1(:,1),E1(:,2)); xlabel('t');ylabel('y');


% Modified Traezodial Method
% param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
% - start_t,end_t : start and end limits of time
% - y0 : Initial condition, y(start_t)=y0
% - stepN : stepN * h = (end_t-start_t) , total EulerExp steps performed
% @ return_MT : Matrix of results. return_MT(:,1) is time t,
% return_MT(:,2) is value y at time t.
function return_MT = ModifiedTraezodial(f,start_t,end_t,y0,stepN)
    % Define variables 
    h = (end_t-start_t)/stepN; %step size
    T = linspace(start_t,end_t,stepN+1); % create&fill T which is step time
    Yt = zeros(1,stepN+1);
    % Start-up condition 
    Yt(1) = y0; 
    % Each loop comes with approximate value of Yt(i+1)
    for i = 1:stepN
        % Yt(i)==Y_i , CurrY==Yt(i+1)==Y_(i+1) 
        % Y_(i+1)= Y_i+h*f(t_i,Y_i) in mathematics. Get estimated value of Y_(i+1)
        CurrY = Yt(i)+h*f(T(i),Yt(i));  
        % Y_(i+1)= Y_i+h*0.5*[f(t_i,Y_i)+f(t_(i+1),Y_(i+1))]
        CurrY = Yt(i)+h*0.5*(f(T(i),Yt(i))+f(T(i+1),CurrY)); % That's approximate value of Y_(i+1)
        Yt(i+1) = CurrY; % Log
    end
    % reture time and Y
    return_MT = [T' Yt'];
end

end