f = @(t,y)y*(2/t+1);
gf = @(t,y)(2/t+1);
E3 = EulerImplicit_Steps(f,gf,1,3,0.3679,400,1e-8,100);
E2 = EulerImplicit_Steps(f,gf,1,3,0.3679,400,1e-18,10); 
% E2 will come up with an error because tolerance is too small.


% Euler Implicit Method (Newton-Raphson Method based) , using total steps. 
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
function return_EIS = EulerImplicit_Steps(f,gf,start_t,end_t,start_value,total_steps,tolerance,maximum_iteration)
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








