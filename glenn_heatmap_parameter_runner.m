function glenn_heatmap_parameter_runner(phiortau,distribution,par2test,discs,start)

% Function that generates coevolution data from numerous serial transfer experiments performed with different parameter settings in the parameter reduced system, using the glenn_epoch_runner_reduced() functions and glenn_epoch_solver_reduced_mevo() functions underneath.

% glenn_epoch_runner_reduced().
    %{
        This function runs a serial transfer experiment for the reduced
        parameter system, over a range of values for one out of three
        parameters. It then outputs a heat-map series construed from the average
        output of the requested epochs and plots mean and modus of the
        phenotype distribution in a line-plot.

        The input it requires:

            phiortau        = first name-flag given to output files
            distribution    = second name-flag given to output files
            par2test        = "K"/"gamma"/"alpha"
                                Determines which of these parameters is
                                varied between the serial transfer
                                experiments.
            discs           = number of serial transfer experiments done
                                (the iteration which corresponds to the
                                maximal parameter value of par2test)
            start           = number of serial transfer experiment to start
                                (the iteration which corresponds to the
                                first parameter value of par2test you wish
                                to create a heatmap for)

        The function refers to 
            glenn_epoch_runner_reduced()
        this function is similar glenn_epoch_runner, but is built to work
        with the reduced parameter set, and the equations that it solves
        are given in glenn_solver_reduced(), which thus follows the reduced
        set of equations.
    %}

    % Initialize standard set:
    birth       = 1;
    death       = .3;
    resurge     = .4;
    infectivity = .25;
    burstsize   = 4;
    decay       = .5;
    fmut        = 1e-3;
    bmut        = 1e-3;
    
    % Make parameter conversions:
    K           = (birth-death)/decay;
    alpha       = resurge/decay;
    beta        = (decay*decay)/birth;
    gamma       = burstsize*infectivity/birth;
    fm          = fmut*infectivity/decay;
    bm          = bmut*infectivity/decay;
    
    % Make parameter ranges for testing:
    gammaset    = linspace((K-alpha),K,discs);
    alphaset    = linspace((K-gamma),K,discs);
    Kset        = linspace(gamma,(gamma+alpha),discs);
    
    
    % Initialize run:
    philist     = linspace(0,1,10);
    taulist     = linspace(0,0.12,10);
    tolrnc      = 1e-6;
    num         = length(philist)*length(taulist);
    eplen       = 150;
    steps       = 1;
    epochmax    = 1e5;
    
    % Make initial strain distribution
    ins(1)=.7;
    ins(2:(1+num))=0;
    ins((2+num):(1+2*num))=0.001/num;
    
    
    % Actually run:
    for ii=start:discs
        if par2test == "K"
            K=Kset(ii);
        else
            if par2test == "gamma"
                gamma=gammaset(ii);
            else
                if par2test == "alpha"
                    alpha=alphaset(ii);
                end
            end
        end
        
        % Calculate:
        if gamma > K-alpha && gamma < K
            string=" -- IN REGIME: YES";
        else
            string=" -- IN REGIME: NO";
        end
        index=num2str(ii,'%i');
        disp("WORKING ON "+par2test+index+string)
        
        [tf,cf,tn,cn,te,~,avg,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner_reduced(eplen,steps,epochmax,philist,taulist,K,alpha,beta,gamma,fm,bm,ins,tolrnc);
                        
        % Save data:
        dlmwrite(phiortau+"_"+distribution+"_"+par2test+"_"+"run_"+index+"_fc"+".txt",cf);
        dlmwrite(phiortau+"_"+distribution+"_"+par2test+"_"+"run_"+index+"_ft"+".txt",tf);
        dlmwrite(phiortau+"_"+distribution+"_"+par2test+"_"+"run_"+index+"_nc"+".txt",cn);
        dlmwrite(phiortau+"_"+distribution+"_"+par2test+"_"+"run_"+index+"_nt"+".txt",tn);
        dlmwrite(phiortau+"_"+distribution+"_"+par2test+"_"+"run_"+index+"_et"+".txt",te);
        dlmwrite(phiortau+"_"+distribution+"_"+par2test+"_"+"run_"+index+"_av"+".txt",avg);

    end
    folder="";
    glenn_heatmap_prepper(philist,taulist,folder,phiortau,distribution,par2test,1:discs)
    glenn_alternate_mevo_plotter(folder,phiortau,distribution,par2test,1:discs,philist,taulist)
end