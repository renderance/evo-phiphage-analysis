function x = glenn_solver_mevo(br,dr,al,be,ga,de,phis,taus,mf,mt,time,x0)
% Single epoch trajectory solver allowing for co-evolution of threshold and lysogeny propensity.
    %{
     glenn_solver_mevo() 

     This function solves competition experiments between temperate bacteriophage 
     of both lysogeny regulating and non-regulating kind. It creates phenotypes with
     combined features from the lists that are provided to it, where one gives the
     maximal propensity for lysogeny and the other provides the regulation threshold 
     of the phenotype. For the integration of the system of ordinary differential 
     equations, it uses either the built-in solver ode45, or ode23s if that 
     fails (though it shouldn't).

       Usage: glenn_solver(
               host birthrate, 
               host deathrate, 
               resurgence rate, 
               burstsize, 
               phage decayrate,  
               list of maximal lysogeny propensities, 
               list of regulator tresholds, 
               forward mutation rate, 
               backward mutation rate, 
               ODE45 list of intervals of integration, 
               starting conditions
           )

     Note that this solver was built for the co-evolution of both the
     maximal propensity for lysogeny and the regulation threshold, in the
     full parameter system. It creates a vector for plist*tlist
     phenotypes.

    %}


    nf = length(phis);
    nt = length(taus);

    % Check if the starting conditions vector has the correct length.
    if size(x0,2)~= 1+nf*nt*2
        fprintf(['Error in glenn_solver initiation: initial vector should have length ',num2str(1+nf*nt*2,'%i'),' according to the phi and tau lists provided.\n'])
        return
    end

    % Solve the system and print which solver was used, in case one failed.
    lastwarn('');
    [~,x] = ode45(@RHS,time,x0);
    [~,msgid]=lastwarn;
    if msgid ~= 0
        [~,x] = ode23s(@RHS,time,x0);
        fprintf(' finished with ODE23S.\n')
    else
        fprintf(' finished with ODE45.\n')
    end

    % Build the input vector of phenotypes.
    function fx = RHS(~,x)

        % Go over the input, set all negative values to zero.
        %{
            An ugly fix to make sure the sytem does not explode in rare cases
            of regulator miscalculation. Should not affect results much
            otherwise.
        %}
        x(x<0)=0;
        
        fx = zeros(1+2*nf*nt,1);

        % The index for sensitive population is one.
        inS = 1;

        % Build the function for the sensitives.
        fx(inS) = br*x(1)*(1-x(1)-sum(x(2:(1+nf*nt))))-dr*x(1)-be*x(1)*sum(x((2+nf*nt):(1+2*nf*nt)));        
        
        % Begin building lysogens.
        for ii = 1:nt
            for jj = 1:nf
                
                % Calculate the current lysogeny propensity of the
                % bacteriophage under examination.
                PHI = heaviside(be*x(1)*sum(x((2+nf*nt):(1+2*nf*nt)))-taus(ii))*phis(jj);
                
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
                
                % Check possibility of mutation, such that corners of our
                % phenotype matrix do not mutate to mutant phenotypes that
                % lie outside of it.
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
                
                
                % Making lysogen function
                fx(inJ) = br*x(inJ)*(1-x(1)-sum(x(2:(1+nf*nt))))-dr*x(inJ)-al*x(inJ)+PHI*be*x(1)*x(inP);

                % Making phage function:
                
                % Add growth of self but substract loss through mutation.
                fx(inP) = (1-omff-omft-ombf-ombt)*ga*(al*x(inJ)+(1-PHI)*be*x(1)*x(inP));
                
                % Add acquisition through forward mutation from lower tau.
                if ii ~= 1
                    fx(inP) = fx(inP)+imft*ga*(al*x(inJ-nf)+(1-PHI)*be*x(1)*x(inP-nf));
                end
                
                % Add acquisition through backward mutation from higher tau.
                if ii ~= nt
                    fx(inP) = fx(inP)+imbt*ga*(al*x(inJ+nf)+(1-PHI)*be*x(1)*x(inP+nf));
                end
                
                % Add acquisition through forward mutation from lower phi.
                if jj ~= 1
                    fx(inP) = fx(inP)+imff*ga*(al*x(inJ-1.)+(1-PHI)*be*x(1)*x(inP-1.));
                end
                
                % Add acquisition through backward mutation from higher phi.
                if jj ~= nf
                    fx(inP) = fx(inP)+imbf*ga*(al*x(inJ+1.)+(1-PHI)*be*x(1)*x(inP+1.));
                end

                % Subtract decay of phage.
                fx(inP) = fx(inP)-de*x(inP);
                
            end
        end
    end
end
                    
                                    
