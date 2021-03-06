function [matrix,xlist,ylist,xtext,ytext,signal] = glenn_heatmap_runner(phiortau,distribution,matrixname,xlistname,ylistname,signalname,identifier,eplist)

% Function that generates heatmap data from numerous serial transfer experiments performed with various epoch lenghts, using the glenn_epoch_runne*() functions and glenn_epoch_solve*() functions underneath.


%{
    This function is magical. No, seriously, it is.

    It performs a number of serial transfers for the evolution of either
    regulators, non-regulators or the co-evolotion of the regulation
    threshold and the maximal lysogeny propensity. It does so for 50
    phenotypes, or in the case of the co-evolution, for 100 (10x10) phenotypes.
    
    It does so with the standard parameter set and performs the serial
    transfer experiment for a range of (mean) epoch lengths.

    Each epoch is calculated by an epoch runner function, which prepare input
    for the solver functions. The output of the serial transfer experiments
    are written to files.

    These files can be read and treated with glenn_heatmap_prepper().

    Usage:
        glenn_heatmap_transfer_runner(
                phiortau        =   "phi"   non-regulators are evolved with
                                            lysogeny propensities from 0 to 1.
                                    "tau"   regulators are evolved with maximal
                                            lysogeny propensity 1 and regulation
                                            thresholds ranging from 0 to 0.11.
                                    "duo"   phenotypes are created with maximal
                                            lysogeny propensities ranging from 0 to
                                            1 and with regulation thresholds
                                            ranging from 0 to 0.12.
                distribution    =   "norm"  the epochs have a constant value
                                    "expd"  the epoch length is drawn from
                                            an exponential distribution.
                                    "mevo"  the epochs have a constant
                                            value but the morevo solver is
                                            used to allow mutation in two
                                            directions.
                matrixname      =   ?       gives the file name flag to which
                                            the eventual heatmap matrix should 
                                            be saved.
                xlistname       =   ?       gives the file name flag to which
                                            the x-axis values should be saved.
                ylistname       =   ?       idem, for the y-axis.
                signalname      =   ?       gives the file name flag to which
                                            the signal should be saved.
                identifier      =   ?       gives the file name flag for
                                            the entire run.
                    SUCH THAT THE  FILE NAMES BECOME:
                        phiortau_distribution_identifier_objectflag.txt

                eplist          =   [...]   a vector containing the epoch lengths
                                            that a serial transfer experiment
                                            should be performed for.
            )

    Note that this script refers to four different types of epoch runners,
    depending on which of the combinations were chosen. It refers to:

        glenn_epoch_runner.m
        glenn_epoch_runner_exprnd.m     draws eplen from exponential distribution
        glenn_epoch_runner_transfer.m   draws transfer size from a lognormal distribution
        glenn_epoch_runner_morevo.m     mutates in two evolvable parameter space directions

%}

    % Initialize parameters:
    birth = 1;
    death = .3;
    resurge = .4;
    infectivity = .25;
    burstsize = 4;
    decay = .5;
    fomut = 1e-3;
    bamut = 1e-3;
    steps = 1;
    maxeps = 1e5;
    
    % Initialize run:
    if phiortau=="phi"
        matrix=[];
        signal=[];
        philist=linspace(0,1,50);
        ins(1)=1-death/birth;
        ins(2:(1+length(philist)))=0;
        ins((2+length(philist)):(1+2*length(philist)))=0.01/50;
        taulist=[];
        xtext="Normalized distribution of \phi's in population";
        ytext="length of one epoch";
        tolrnc=1e-7;
    else
        if phiortau=="tau"
            matrix=[];
            signal=[];
            philist=[];
            taulist=linspace(0,.11,50);
            ins(1)=1-death/birth;
            ins(2:(1+length(taulist)))=0;
            ins((2+length(taulist)):(1+2*length(taulist)))=0.01/50;
            xtext="Normalized distribution of \tau's in population";
            ytext="length of one epoch";
            tolrnc=1e-6;
        else
            if phiortau=="duo"
                matrix=[];
                signal=[];
                philist=linspace(0,1,10);
                taulist=linspace(0,0.12,10);
                ins(1)=1-death/birth;
                ins(2:(1+length(philist)*length(taulist)))=0;
                ins((2+length(philist)*length(taulist)):(1+2*(length(philist)*length(taulist))))=0.01/(length(philist)*length(taulist));
                xtext="Normalized distribution of \tau's or \phi's in population";
                ytext="length of one epoch";
                tolrnc=1e-6;
            end
        end
    end
    
    % Actually run:
    for eplen = eplist
        if eplen>10000
            maxeps=1;
        end
        
        % Calculate:
        fprintf(['\nWORKING ON EPOCH LENGTH ',num2str(eplen,'%i'),'\n']);
        if distribution == "norm"
            [tf,cf,tn,cn,te,ce,avg,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner(eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,ins,tolrnc);
        else
            if distribution == "expd"
                [tf,cf,tn,cn,te,ce,avg,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner_exprnd(eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,ins,tolrnc);
            else
                if distribution == "tran"
                    [tf,cf,tn,cn,te,ce,avg,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner_transfer(eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,ins,tolrnc);
                else
                    if distribution == "mevo"
                        [tf,cf,tn,cn,te,ce,avg,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner_morevo(eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,ins,tolrnc);
                        xtext="Distribution of \tau's in population";
                        ytext="Distribution of \phi_{max} in population";
                    end
                end
            end
        end
        
        % Save data:
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+"run_"+num2str(eplen,'%i')+"_fc"+".txt",cf);
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+"run_"+num2str(eplen,'%i')+"_ft"+".txt",tf);
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+"run_"+num2str(eplen,'%i')+"_nc"+".txt",cn);
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+"run_"+num2str(eplen,'%i')+"_nt"+".txt",tn);
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+"run_"+num2str(eplen,'%i')+"_et"+".txt",te);
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+"run_"+num2str(eplen,'%i')+"_av"+".txt",avg);
        
        if distribution ~= "mevo"
            if eplen>10000
                matrix=cat(1,matrix,mean(cn(1+length(philist)+length(taulist)+1:end,end)+cn(2:1+length(philist)+length(taulist))));
            else
                avgn=zeros(length(1+length(philist)+length(taulist)+1:1+2*(length(philist)+length(taulist))),1);
                for uu=1:(length(philist)+length(taulist))
                    ux=(1+length(philist)+length(taulist)+uu);
                    avgn(uu,1)=mean(cn(ux,(end-1000*eplen*1/steps):end)+cn(1+uu,(end-1000*eplen*1/steps):end));
                end
                matrix=cat(1,matrix,avgn(1:end,1)');
                if isempty(taulist)==0
                    signal=cat(1,signal,infectivity*sum(ce(1+length(philist)+length(taulist)+1:end,end)*ce(1,end)));
                end
            end
        else
            avgn=zeros(length(1+length(philist)*length(taulist)+1:1+2*(length(philist)*length(taulist))),1);
            for uu=1:(length(philist)*length(taulist))
                ux=(1+length(philist)*length(taulist)+uu);
                avgn(uu,1)=mean(cn(ux,(end-1000*eplen*1/steps):end)+cn(1+uu,(end-1000*eplen*1/steps):end));
            end
            avgn2=zeros(length(philist),length(taulist));
            for vv = 1:length(philist)
                for ww = 1:length(taulist)
                    avgn2(vv,ww,1)=avgn((ww-1)*length(philist)+vv,1);
                end
            end
            matrix=cat(3,matrix,avgn2(:,:,1));
        end
    end
    
    % Output data files:
    if distribution ~= "mevo"
        % Store axes tick values:
        if (phiortau=="phi")
            xlist=philist;
            ylist=eplist;
        else
            if phiortau=="tau"
                xlist=taulist;
                ylist=eplist;
            else
                if phiortau=="duo"
                    xlist=[philist taulist];
                    ylist=eplist;
                end
            end
        end
        % Store axes ticks:
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+xlistname+".txt",xlist);
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+ylistname+".txt",ylist);
        % Store matrix of data:
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+matrixname+".txt",matrix);
        % Store matrix of signal:
        if isempty(signal)==0
            dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+signalname+".txt",signal);
        end
    else
        % Store axes tick values:
        xlist=taulist;
        ylist=philist;
        % Store axes ticks:
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+xlistname+".txt",xlist);
        dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+ylistname+".txt",ylist);
        % Store matrixes of data:
        for uu = 1:length(eplist)
            dlmwrite(phiortau+"_"+distribution+"_"+identifier+"_"+matrixname+"_"+num2str(eplist(uu),'%i')+".txt",matrix(:,:,uu));
        end
    end
end