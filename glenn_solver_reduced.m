function x = glenn_solver_reduced(K,alpha,gamma,rho,plist,tlist,fmut,bmut,time,x0)

% Solver for the trajectory of one epoch of the parameter reduced system

    %{
     glenn_solver_reduced()

     This function solves competition experiments between temperate bacteriophage 
     of both lysogeny regulating and non-regulating kind. It determines which 
     to add of either from the lists that are provided to it. For the integration
     of the system of ordinary differential equations, it uses either the
     built-in solver ode45, or ode23s if that fails.

       Usage: glenn_solver(
               maximal host organism density parameter, 
               adjusted resurgence rate parameter, 
               adjusted burstsize parameter, 
               signal scaling factor,
               list of non-regulator lysogeny fractions, 
               list of regulator tresholds, 
               forward mutation rate, 
               backward mutation rate, 
               ODE45 interval of integration, 
               starting conditions
           )
    %}

% Check if the starting vector is as long as it should be:

num=size(list,2);
if size(x0,2)~= 1+num*2
    fprintf(['Error in glenn_solver initiation: initial vector should have length ',num2str(1+num*2,'%i'),' according to the phi and tau lists provided.\n'])
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
        
        phis = zeros(plist,tlist);
        phis = heaviside(rho*x(1)*sum(x((1+num+1):(1+num+num)))-tlist(:)).*plist(:);

            
        
        % Create an empty column-vector for all the equations:
        fx = zeros(1+2*num,1);
        
        % Create the first index, which always contains the sensitives:
        fx(1) =                         x(1).*(K-x(1)-sum(x(2:(1+num))))...
                                        -x(1).*sum(x((1+num+1):end));
                
        % Create the lysogens:
        fx(2:(1+num)) =                x(2:(1+num)).*(K-x(1)-sum(x(2:(1+num)))-alpha)...
                                        +phis(:).*x(1).*x((1+num+1):end);
                                    
        % Create the first phage:
        fx(1+num+1) =                  gamma*(alpha*x(2)+(1-phis(1))*x(1)*x(1+num+1))*(1-fmut)+...
                                        gamma*(alpha*x(3)+(1-phis(2))*x(1)*x(1+num+2))*bmut...
                                        -x(1+num+1);
                                    
        % Create the last phage:
        fx(1+num+num) =               gamma*(alpha*x(1+num)+(1-phis(end))*x(1)*x(1+num+num))*(1-bmut)+...
                                        gamma*(alpha*x(num)+(1-phis(end-1))*x(1)*x(num+num))*fmut...
                                        -x(1+num+num);
                                    
        % Create the phages that can mutate in both directions:
        fx((1+num+2):(num+num)) =    gamma*(alpha*x((3):(num))+(1-phis(2:(end-1)))*x(1)*x((1+num+2):(num+num)))*(1-bmut-fmut)+...
                                        gamma*(alpha*x(((3):(num))-1)+(1-phis((2:(end-1))-1))*x(1)*x(((1+num+2):(num+num))-1))*fmut+...
                                        gamma*(alpha*x(((3):(num))+1)+(1-phis((2:(end-1))+1))*x(1)*x(((1+num+2):(num+num))+1))*bmut...
                                        -x((1+num+2):(num+num));
                                    
    end
end