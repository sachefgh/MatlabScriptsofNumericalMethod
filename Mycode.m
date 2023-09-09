f=@(t,a)[0,1;-2,-0.1]*a+[0;0.2*cos(2*t)+0.5*sin(t)];

interval= 0.1;
E1 = ModifiedTrapezodial_MAT_Interval(f,0,50,[0.1;0],interval);
E2 = Runge_Kutta_MAT_Interval(f,0,50,[0.1;0],interval);
%xstring = sprintf('t (h=%g)',interval);
plot(E1(1,:),E1(2,:),'--','Color','red'); hold;
%xlabel(xstring,'Color','magenta');ylabel('Y','Color','magenta');
plot(E2(1,:),E2(2,:),'--','Color','blue');

% Modified Trapezodial Method - MATRIX - using total steps
% This is to solve 2nd order IVP INDIRECTLY, by modify all relations
% between multiple order into function f.
% param:
% - f : refers to function handle of f(t,[a]) |||  a/dt=f(t,[a]) 
%       @@ note that [a] should be a n*1 array
%       f = [a]/dt = [K][a] + [F] ;[K] is constant matrix, [F] is function
%       matrix related to time t
%     e.g.  f=@(t,a)[0,1;-2,-0.1]*a+[0;0.2*cos(2*t)+0.5*sin(t)];
% - start_t,end_t : start and end limits of time
% - a0 : Initial condition, a(start_t)= [a0]
%       @@ note that [a0] should be a n*1 array
% - stepN : stepN * h = (end_t-start_t) , total EulerExp steps performed
%
% @ return_MTM : Matrix of results. First row return_MTM(1,i) is time t_i,
%   return_MTM(2:end,i)-->is [a]_i (n*1 array) at time t_i.
%
function return_MTM = ModifiedTrapezodial_MAT_steps(f,start_t,end_t,a0,stepN)
    % define variables
    T = linspace(start_t,end_t,stepN+1); 
    a = cell(1, stepN+1); % using cell array to store value a; a{1,i+1}->[a]_i 
    h = (end_t-start_t)/stepN; %step size
    % Init condition
    a{1,1} = a0;
    % Iteration loop, using Modified Traezodial Method
    for i = 1:stepN
        Curr_a = a{1,i} + h*f(T(i),a{1,i}); % Get estimated value of [a_(i+1)]
        Curr_a = a{1,i} + 0.5*h*(f(T(i),a{1,i})+f(T(i+1),Curr_a));% That's approximate value of [a_(i+1)]
        a{1,i+1} = Curr_a; % store value [a] for current iteration
    end
    % return value
    % IMPORTANT : should covert cell a 2 Mat for proper output!!
    return_MTM = T;  % fill first row with time t
    return_MTM = [return_MTM;cell2mat(a)]; % merge vertically
end

% Modified Trapezodial Method - MATRIX - using step Interval
% This is to solve 2nd order IVP INDIRECTLY, by modify all relations
% between multiple order into function f.
% param:
% - f : refers to function handle of f(t,[a]) |||  a/dt=f(t,[a]) 
%       @@ note that [a] should be a n*1 array
%       f = [a]/dt = [K][a] + [F] ;[K] is constant matrix, [F] is function
%       matrix related to time t
%     e.g.  f=@(t,a)[0,1;-2,-0.1]*a+[0;0.2*cos(2*t)+0.5*sin(t)];
% - start_t,end_t : start and end limits of time
% - a0 : Initial condition, a(start_t)= [a0]
%       @@ note that [a0] should be a n*1 array
% - h : stepN * h = (end_t-start_t) , interval between tˇi+1 and tˇi
%
% @ return_MTM : Matrix of results. First row return_MTM(1,i) is time t_i,
%   return_MTM(2:end,i)-->is [a]_i (n*1 array) at time t_i.
%
function return_MTM = ModifiedTrapezodial_MAT_Interval(f,start_t,end_t,a0,h)
    stepN = ceil((end_t-start_t)/h);
    % define variables
    T = zeros(1,stepN+1);
    a = cell(1, stepN+1);
    % fill time array.
    T(1)=start_t; T(stepN+1)=end_t;
    for j=1:(stepN-1)
        T(j+1)=T(j)+h;
    end
    % Init condition
    a{1,1} = a0;    
    % Iteration loop, using Modified Traezodial Method
    for i = 1:stepN
        Curr_a = a{1,i} + h*f(T(i),a{1,i}); % Get estimated value of [a_(i+1)]
        Curr_a = a{1,i} + 0.5*h*(f(T(i),a{1,i})+f(T(i+1),Curr_a));% That's approximate value of [a_(i+1)]
        a{1,i+1} = Curr_a; % store value [a] for current iteration
    end
    % return value
    % IMPORTANT : should covert cell a 2 Mat for proper output!!
    return_MTM = T;  % fill first row with time t
    return_MTM = [return_MTM;cell2mat(a)]; % merge vertically

end

% Runge-Kutta Forth-order(RK4) Method - non-matrix - using total steps 
% param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
% - start_t,end_t : start and end limits of time
% - y0 : Initial condition, y(start_t)=y0
% - stepN : stepN * h = (end_t-start_t) , total EulerExp steps performed
% @ return_RK : Matrix of results. return_MT(:,1) is time t,
% return_RK(:,2) is value y at time
function return_RK = Runge_Kutta_steps(f,start_t,end_t,y0,stepN)
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

% Runge-Kutta Forth-order(RK4) Method - MATRIX - using total steps
% This is to solve 2nd order IVP INDIRECTLY, by modify all relations
% between multiple order into function f.
% param:
% - f : refers to function handle of f(t,[a]) |||  a/dt=f(t,[a]) 
%       @@ note that [a] should be a n*1 array
%       f = [a]/dt = [K][a] + [F] ;[K] is constant matrix, [F] is function
%       matrix related to time t. 
%     e.g.  f=@(t,a)[0,1;-2,-0.1]*a+[0;0.2*cos(2*t)+0.5*sin(t)];
% - start_t,end_t : start and end limits of time
% - a0 : Initial condition, a(start_t)= [a0]
%       @@ note that [a0] should be a n*1 array
% - stepN : stepN * h = (end_t-start_t) , total  steps performed
% @ return_RK : Matrix of results. First row return_RK(1,i) is time t_i,
%   return_RK(2:end,i)-->is [a]_i (n*1 array) at time t_i.
%
function return_RK = Runge_Kutta_MAT_steps(f,start_t,end_t,a0,stepN)
    % define variables
    T = linspace(start_t,end_t,stepN+1); 
    a = cell(1, stepN+1); % using cell array to store value a; a{1,i+1}->[a]_i 
    h = (end_t-start_t)/stepN; %step size
    % Init condition
    a{1,1} = a0;
    % Iteration using RK4 Method
    for i=1:stepN
        ai05_1 = a{1,i}+0.5*h*f(T(i),a{1,i});
        ai05_2 = a{1,i}+0.5*h*f(T(i)+0.5*h,ai05_1);
        ai1_3 = a{1,i}+h*f(T(i)+0.5*h,ai05_2);
        a{1,i+1} =  a{1,i} + (1/6)*h*(  f(T(i),a{1,i}) + 2*f(T(i)+0.5*h,ai05_1) + 2*f(T(i)+0.5*h,ai05_2) + f(T(i+1),ai1_3) );
    end
    % return value
    % IMPORTANT : should covert cell a 2 Mat for proper output!!
    return_RK = T; % fill first row with time t
    return_RK = [return_RK;cell2mat(a)]; % merge vertically
end

% Runge-Kutta Forth-order(RK4) Method - MATRIX - using Interval
% This is to solve 2nd order IVP INDIRECTLY, by modify all relations
% between multiple order into function f.
% param:
% - f : refers to function handle of f(t,[a]) |||  a/dt=f(t,[a]) 
%       @@ note that [a] should be a n*1 array
%       f = [a]/dt = [K][a] + [F] ;[K] is constant matrix, [F] is function
%       matrix related to time t. 
%     e.g.  f=@(t,a)[0,1;-2,-0.1]*a+[0;0.2*cos(2*t)+0.5*sin(t)];
% - start_t,end_t : start and end limits of time
% - a0 : Initial condition, a(start_t)= [a0]
%       @@ note that [a0] should be a n*1 array
% - h : stepN * h = (end_t-start_t) , interval
% @ return_RK : Matrix of results. First row return_RK(1,i) is time t_i,
%   return_RK(2:end,i)-->is [a]_i (n*1 array) at time t_i.
%
function return_RK = Runge_Kutta_MAT_Interval(f,start_t,end_t,a0,h)
    stepN = ceil((end_t-start_t)/h);
    % define variables
    T = zeros(1,stepN+1);
    a = cell(1, stepN+1);
    % fill time array.
    T(1)=start_t; T(stepN+1)=end_t;
    for j=1:(stepN-1)
        T(j+1)=T(j)+h;
    end
     % Init condition
    a{1,1} = a0;
    % Iteration using RK4 Method
    for i=1:stepN
        ai05_1 = a{1,i}+0.5*h*f(T(i),a{1,i});
        ai05_2 = a{1,i}+0.5*h*f(T(i)+0.5*h,ai05_1);
        ai1_3 = a{1,i}+h*f(T(i)+0.5*h,ai05_2);
        a{1,i+1} =  a{1,i} + (1/6)*h*(  f(T(i),a{1,i}) + 2*f(T(i)+0.5*h,ai05_1) + 2*f(T(i)+0.5*h,ai05_2) + f(T(i+1),ai1_3) );
    end
    % return value
    % IMPORTANT : should covert cell a 2 Mat for proper output!!
    return_RK = T; % fill first row with time t
    return_RK = [return_RK;cell2mat(a)]; % merge vertically
end

% Modified Traezodial Method - non-matrix - using total steps
% param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
% - start_t,end_t : start and end limits of time
% - y0 : Initial condition, y(start_t)=y0
% - stepN : stepN * h = (end_t-start_t) , total EulerExp steps performed
% @ return_MT : Matrix of results. return_MT(:,1) is time t,
% return_MT(:,2) is value y at time t.
function return_MT = ModifiedTraezodial_steps(f,start_t,end_t,y0,stepN)
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

% Euler Implicit Method (Newton-Raphson Method based) - non-matrix - using total steps. 
% Param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
% - gf : refers to function handle of [əf(t,y)/əy] ; 
%           e.g. -------->    f=@(t,y)(2/t+1)
%
% - start_t,end_t : start and end limits of time
% - start_value : Initial condition, y(start_t)=start_value
% - total_steps : total_steps*h = (end_t-start_t)
% - tolerance : tolerance of [yˇi+1-yˇi-hf(tˇi+1,yˇi+1)].  e.g. ----> 1e-8
% - maximum number of iterations before break&throw an error.
% @ return_EIS : return Matrix of results. return_EIS(:,1) is time t,
% return_EIS(:,2) is value y at time t.
%////////////////////////////////////////////////////////////////////%
% // There's a document attached to explain why I think by providing//
% gf and f can function well, from my deduction.                    //
%////////////////////////////////////////////////////////////////////%
function return_EIS = EulerImplicit_steps(f,gf,start_t,end_t,start_value,total_steps,tolerance,maximum_iteration)
    h = (end_t-start_t)/total_steps;
    T = linspace(start_t,end_t,total_steps+1);
    Yt = zeros(1,total_steps+1);
    Yt(1) = start_value;
    for m=1:total_steps 
        PrevY = Yt(m); % Read Yi
        CurrY = PrevY;    % Set Yˇi+1 = Yˇi (known) at the beginning.
        erro_flag = true;
        for n= 1:maximum_iteration
        % If offset smaller tolerance occur, make iteration once at least.
            CurrY = CurrY - (CurrY-PrevY-h*f(T(m+1),CurrY))/(1-h*gf(T(m+1),CurrY)); % Delta and F(Yˇi+1) refreshed.
            if abs(CurrY-PrevY-h*f(T(m+1),CurrY)) < tolerance % exit loop while Yˇi+1  converge.
                erro_flag = false;
                break
            end
        end
        % Now check if operation successful.
        if erro_flag 
            error('Exceeding maximum number of iterations!! FORCE STOP');
        end 
        Yt(m+1) = CurrY; % store Yˇi+1 
    end
    return_EIS = [T' Yt'];
end

% Euler Explicit Method - non-matrix - using total steps.
% Param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
%
% - start_t,end_t : start and end limits of time
% - start_value : Initial condition, y(start_t)=start_value
% - total_steps : total_steps*h = (end_t-start_t)
% @ return_EES : return Matrix of results. return_EES(:,1) is time t,
% return_EES(:,2) is value y at time t.
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

% Euler Explicit Method - non-matrix - using step interval.
% Param:
% - f : refers to function handle of f(t,y) |||  y'=f(t,y) ||| yˇi+1=yˇi + hf(tˇi,yˇi)
%           e.g. -------->    f=@(t,y)y*(2/t+1);
%
% - start_t,end_t : start and end limits of time
% - start_value : Initial condition, y(start_t)=start_value
% - h : interval between tˇi+1 and tˇi
% @ return_EEI : return Matrix of results. return_EEI(:,1) is time t,
% return_EEI(:,2) is value y at time t.
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

