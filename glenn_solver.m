function x = glenn_solver(b,d,alpha,beta,gamma,delta,phis,taus,fmut,bmut,time,x0)
% Single epoch trajectory solver for full system for regulators and non-regulators.
    %{
     glenn_solver() 
     This function solves competition experiments between temperate bacteriophage 
     of both lysogeny regulating and non-regulating kind. It determines which 
     to add of either from the lists that are provided to it. For the integration
     of the system of ordinary differential equations, it uses either the
     built-in solver ode45, or ode23s if that fails.

       Usage: glenn_solver(
               host birthrate, 
               host deathrate, 
               resurgence rate, 
               burstsize, 
               phage decayrate, 
               list of non-regulator lysogeny fractions, 
               list of regulator tresholds, 
               forward mutation rate, 
               backward mutation rate, 
               ODE45 interval of integration, 
               starting conditions
           )
    %}

% Check if the starting vector is as long as it should be:

numf=size(phis,2);
numt=size(taus,2);
if size(x0,2)~= 1+numf*2+numt*2
    fprintf(['Error in glenn_solver initiation: initial vector should have length ',num2str(1+numf*2+numt*2,'%i'),' according to the phi and tau lists provided.\n'])
    return
end


% Invoke the ode45 solver for the system of equations, and if it fails, try ode23s:
lastwarn('');
[t,x] = ode45(@RHS,time,x0);
[~,msgid]=lastwarn;
    if msgid ~= 0
        [t,x] = ode23s(@RHS,time,x0);
        fprintf(' finished with ODE23S.\n')
    else
        fprintf(' finished with ODE45.\n')
    end
    
    
    
    % Construct the Right-Hand-Side vector of equations used by the solver:
    function fx = RHS(t,x)
        
        x(x<0)=0;
        
        % Create an empty column-vector for all the equations:
        fx = zeros(1+2*numf+2*numt,1);
        
        % Create the first index, which always contains the sensitives:
        fx(1) =                         b*x(1)*(1-sum(x(1:(numf+numt+1))))...
                                        -d*x(1)...
                                        -beta*x(1)*sum(x((1+numf+numt+1):(1+2*numf+2*numt)));
                        
        % Create the lysogens for the non-regulator:
        if numf ~=0
            for num=1:numf
                fx(1+num) =             b*(1-sum(x(1:(numf+numt+1))))*x(1+num)...
                                        +beta*x(1)*phis(num)*x(1+numf+numt+num)...
                                        -d*x(1+num)...
                                        -alpha*x(1+num);
            end
        end
        
        % Create the regulator phi value:
        regulation(:)=heaviside(beta*x(1)*sum(x((1+numf+numt):(1+2*numf+2*numt)))-taus(:));
        
        % Create the lysogens for the regulator:
        if numt ~=0
            for num=1:numt
                fx(1+numf+num) =        b*(1-sum(x(1:(numf+numt+1))))*x(1+numf+num)...
                                        +beta*x(1)*regulation(num)*x(1+numf+numt+numf+num)...
                                        -d*x(1+numf+num)...
                                        -alpha*x(1+numf+num);
            end
        end
        
        % Create the phage growths for the non-regulator, without decay:
        if numf ~=0
            for num=1:numf
                fx(1+numf+numt+num)=    (gamma*beta*x(1)*(1-phis(num))*x(1+numf+numt+num)...
                                        +gamma*alpha*x(1+num));
            end
        end
        
        % Create the phage growths for the regulator, without decay:
        if numt ~=0
            for num=1:numt
                fx(1+2*numf+numt+num)=  (gamma*beta*x(1)*(1-regulation(num))*x(1+2*numf+numt+num)...
                                        +gamma*alpha*x(1+numf+num));
            end
        end
        
        % Make mutation matrix for the non-regulators and multiply:
        if numf ~=0
            Mphi=zeros(numf,numf);
            for ding = 1:numf
                wut = ding-1;
                if wut > 0
                    Mphi(ding,wut)=bmut;
                end
                wut = ding+1;
                if wut < numf+1
                    Mphi(ding,wut)=fmut;
                end
                Mphi(ding,ding)=1-sum(Mphi(ding,:));
            end
        end
        
        % Make mutation matrix for the regulators:
        if numt ~=0
            Mtau=zeros(numt,numt);
            for ding = 1:numt
                wut = ding-1;
                if wut > 0
                    Mtau(ding,wut)=bmut;
                end
                wut = ding+1;
                if wut < numt+1
                    Mtau(ding,wut)=fmut;
                end
                Mtau(ding,ding)=1-sum(Mtau(ding,:));
            end
        end
        
        % Append phage decays to the non-regulator phages, whilst mutating them:
        if numf ~=0
            fx((1+numf+numt+1):(1+numf+numt+numf))= Mphi*fx((1+numf+numt+1):(1+numf+numt+numf))...
                                                   -delta*x((1+numf+numt+1):(1+numf+numt+numf));
        end
        
        % Append phage decays to the regulator phages, agains whilst mutating them:
        if numt ~=0
            fx((1+numf*2+numt+1):(1+numf*2+numt*2))= Mtau*fx((1+numf*2+numt+1):(1+numf*2+numt*2))...
                                                    -delta*x((1+numf*2+numt+1):(1+numf*2+numt*2));
        end
    end


end