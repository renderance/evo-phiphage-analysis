function x = glenn_solver_reduced_mevo(K,alpha,gamma,rho,plist,tlist,mf,mt,time,x0)

% Epoch trajectory solver for parameter reduced system, allowing for co-evolution of lysogeny propensity and threshold.

    %{
     glenn_solver_reduced_mevo() 

     This function solves competition experiments between temperate bacteriophage 
     of both lysogeny regulating and non-regulating kind. It creates phenotypes with
     combined features from the lists that are provided to it, where one gives the
     maximal propensity for lysogeny and the other provides the regulation threshold 
     of the phenotype. For the integration of the system of ordinary differential 
     equations, it uses either the built-in solver ode45, or ode23s if that 
     fails (though it shouldn't).

       Usage: glenn_solver(
               maximal host organism density parameter, 
               adjusted resurgence rate parameter, 
               adjusted burstsize parameter, 
               signal scaling factor, 
               list of maximal lysogeny propensities, 
               list of regulator tresholds, 
               forward mutation rate, 
               backward mutation rate, 
               ODE45 list of intervals of integration, 
               starting conditions
           )

     Note that this solver was built for the co-evolution of both the
     maximal propensity for lysogeny and the regulation threshold, in the
     minimal parameter system. It creates a vector for plist*tlist
     phenotypes.

    %}

    
% Check if the starting vector is as long as it should be:

nf=length(plist);
nt=length(tlist);
num=size(plist,2)*size(tlist,2);
if size(x0,2)~= 1+num*2
    fprintf(['Error in glenn_solver initiation: initial vector should have length ',num2str(1+num*2,'%i'),' according to the phi and tau lists provided.\n'])
    return
end


% Invoke the ode45 solver for the system of equations, and if it fails, try ode23s:
lastwarn('');
[~,x] = ode45(@RHS,time,x0);
[~,msgid]=lastwarn;
    if msgid ~= 0
        [~,x] = ode23s(@RHS,time,x0);
        fprintf(' finished with ODE23S.\n')
    else
        fprintf(' finished with ODE45.\n')
    end

    
    % Construct the Right-Hand-Side vector of equations used by the solver:
    function fx = RHS(~,x)
        
        x(x<0)=0;
        
        % Make a matrix with the right phi values.
        phi = heaviside(rho*x(1)*sum(x((1+num+1):(1+num+num)))-tlist(:)).*plist(:)';
        
        % Create an empty column-vector for all the equations:
        fx = zeros(1+2*num,1);
        
        % Create the first index, which always contains the sensitives:
        fx(1) =                         x(1).*(K-x(1)-sum(x(2:(1+num))))...
                                        -x(1).*sum(x((1+num+1):end));

        for ii = 1:length(tlist)
            for jj = 1:length(plist)
                
                % Lysogen with ii'th tau and jj'th maximal phi index
                inJ = 1+(ii-1)*nf+jj;

                % Corresponding phage index
                inP = inJ+nf*nt;

                % Set mutations for forward and backward
                imft = mt;
                imff = mf;
                imbt = mt;
                imbf = mf;

                omft = mt;
                omff = mf;
                ombt = mt;
                ombf = mf;

                % Check possibility of mutation
                if ii == 1
                    ombt = 0;
                end
                if ii == nt
                    omft = 0;
                end
                if jj == 1
                    ombf = 0;
                end
                if jj == nf
                    omff = 0;
                end
                                    
                % Create the lysogens:
                fx(inJ) = x(inJ)*(K-x(1)-sum(x(2:(1+num)))-alpha)+phi(ii,jj)*x(1)*x(inP);
                                            
                                            
                % Making phage function:
                
                % Add growth of self but substract loss through mutation.
                fx(inP) = (1-omff-omft-ombf-ombt)*gamma*(alpha*x(inJ)+(1-phi(ii,jj))*x(1)*x(inP));
                
                % Add acquisition through forward mutation from lower tau.
                if ii ~= 1
                    fx(inP) = fx(inP)+imft*gamma*(alpha*x(inJ-nf)+(1-phi(ii-1,jj))*x(1)*x(inP-nf));
                end
                
                % Add acquisition through backward mutation from higher tau.
                if ii ~= nt
                    fx(inP) = fx(inP)+imbt*gamma*(alpha*x(inJ+nf)+(1-phi(ii+1,jj))*x(1)*x(inP+nf));
                end
                
                % Add acquisition through forward mutation from lower phi.
                if jj ~= 1
                    fx(inP) = fx(inP)+imff*gamma*(alpha*x(inJ-1)+(1-phi(ii,jj-1))*x(1)*x(inP-1));
                end
                
                % Add acquisition through backward mutation from higher phi.
                if jj ~= nf
                    fx(inP) = fx(inP)+imbf*gamma*(alpha*x(inJ+1)+(1-phi(ii,jj+1))*x(1)*x(inP+1));
                end

                % Subtract decay of phage.
                fx(inP) = fx(inP)-x(inP);
                
            end
        end
    end
end