function [time_full,conc_full,time_norm,conc_norm,time_end,conc_end,avg_ep,eplen,epochstep,epochmax,philist,taulist,final_epoch,K,alpha,beta,gamma,fomut,bamut] = glenn_epoch_runner_reduced(eplen,epochstep,epochmax,philist,taulist,K,alpha,beta,gamma,fmut,bmut,inits,tolrnc)
% Solves simple serial transfer experiment competition equation of phage strains differing in lysogeny rates and decision mechanisms, using glenn_solver_reduced_mevo().


    %%% INITIALIZE %%%
    rng('shuffle')
    runtime=tic;


    %%% PARAMETERS %%%
    tic;
    fprintf('\n\nStarting parameterization...\n')

    % Create host properties

    % Create phage properties
    phis =          philist;
    taus =          taulist;
    num =           length(philist)+length(taulist);

    % Create mutation rates
    fomut = fmut;
    bamut = bmut;

    % Initial strain distribution
    init_x = inits';

    paramtime = toc;
    fprintf(['Time to parameterize system: ',num2str(paramtime,'%.2f'),' seconds.\n'])


    %%% SOLVE THE MODEL %%%
    tic;
    fprintf('Starting solving...\n')

    % Create empty solution objects
    x1 = [];
    x2 = [];
    w = [];
    zz = 0;
    tl1 = [];
    epoch = 0;
    close_nuff = 0;
    dt = epochstep;
    t0 = 0;
    tol = tolrnc;
    vv=0;
    ltstp = 0;
    avgM=[];
    sumM=[];
    te=0;

    % Solve until no longer changes
    while epoch<epochmax && close_nuff < 100*(eplen/dt)*(1+2*num)

        % State which epoch is currently being calculated
        epoch = epoch+1;
        epochlength = eplen;

        % Time vector
        if t0 == 0
            t = t0:dt:epochlength;
            t2 = t(2:end);
            si = 2;
        else
            t = [0 t0:dt:epochlength];
            t2 = t0:dt:epochlength;
            si = 2;
        end
        if t(end) ~= epochlength
            t = [t epochlength];
            ei = length(t)-1;
            t0 = dt-(t(end)-t(end-1));
        else
            t0 = 0;
            ei = length(t);
        end
        te=cat(1,te,te(end)+t(end));

        % Solve the model for an epoch
        fprintf(['Epoch ',num2str(epoch,'%i'),' '])
        y = glenn_solver_reduced_mevo(K,alpha,beta,gamma,phis,taus,fomut,bamut,t,init_x');      
        y(y<0)=0;
        
        % Create a normalized solution matrix for an epoch
        z = zeros(size(y));
        for ii = 1:size(y,1)
            z(ii,1) = 1;
            z(ii,2:(1+num)) = y(ii,2:(1+num))./sum(y(ii,2:(1+num)));
            z(ii,(2+num):end) = y(ii,(2+num):end)./sum(y(ii,(2+num):end));
        end
        z(isnan(z))=0;

        % Concatanate these solutions to the full matrix
        x1 = cat(1,x1,y);
        w = cat(1,w,y(end,:));
        if epoch == 1
            tl1 = t;
            tl2 = t2;
            x2 = z(si:ei,:);
        else
            tl2 = cat(2,tl2,t2+tl1(end));
            tl1 = cat(2,tl1,t+tl1(end));
            x2 = cat(1,x2,z(si:ei,:));
        end
        
        % Create a new starting vector
        init_x(1) = .7;
        init_x(2:(1+num)) = 0;
        init_x((2+num):end) = z(end,(2+num):end) .* sum(inits((2+num):end));
        
        % Create running average matrix
        if epoch > 1000
            sumM = cat(1,sumM,zeros(size(z(si:ei,:))));
            avgM = cat(1,avgM,zeros(size(z(si:ei,:))));
            for uu = (1:length(si:ei))
                vv = ltstp+uu;
                for yy = 1:(1+2*num)
                    if vv == 1
                        sumM(vv,yy)=z((si-1+uu),yy);
                        avgM(vv,yy)=sumM(vv,yy)/vv;
                    else
                        sumM(vv,yy)=sumM(vv-1,yy)+z((si-1+uu),yy);
                        avgM(vv,yy)=sumM(vv,yy)/vv;
                        if abs(avgM(vv,yy)-avgM(vv-1,yy)) < tol %tolrnc*eps(abs(avgM(vv,yy)-avgM(vv-1,yy)))
                            close_nuff = close_nuff+1;
                        else
                            close_nuff = 0;
                        end
                    end
                end
            end
            ltstp = vv;
        end
    end

    solvetime = toc;
    fprintf(['Time to solve system: ',num2str(solvetime,'%.2f'),' seconds.\n'])


    %%% PREPARATION %%%
    tic;
    disp([10,'Starting data preparation...'])

    % Write out the outputs
    time_full=tl1;
    conc_full=x1';
    time_norm=tl2;
    conc_norm=x2';
    time_end=te;
    conc_end=w';
    avg_ep=avgM;
    final_epoch=epoch;

    preptime = toc;
    fprintf(['Time to prepare data: ',num2str(preptime,'%.2f'),' seconds.\n'])

%     figure
%     subplot(1,3,1)
%     plot(avgM(:,1))
%     
%     subplot(1,3,2)
%     plot(avgM(:,2:(1+num)))
%     hold on
%     plot(avgM(:,1),'linestyle',':','color','k')
%     
%     subplot(1,3,3)
%     plot(avgM(:,(2+num):(1+2*num)))
%     hold on
%     plot(avgM(:,1),'linestyle',':','color','k')
    
    %%% FINISH %%%
    runtime=toc(runtime);
    fprintf(['\n\nAll runner operations are completed!\n\n\tCompletion after: ',num2str(runtime,'%.2f'),' seconds.\n\n'])

end