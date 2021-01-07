function glenn_figure_runner(distribution,eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,inits)

% Function that runs a serial transfer experiment with a number of possible settings and uses glenn_figure_awesome() to plot the full trajectory of the experiment as a figure.

%{
    This function runs a serial transfer experiment according to input and
    uses the output generated to create figures. In total, it will generate
    six figures, in two sets of three. The first set works best for
    non-regulators as it does not display the signal. The second set works
    best for regulator phenotypes and combined competition as it does
    contain a panel for the signal.

    The function relies on all four epoch runners, and uses the one defined
    by the first argument.

    glenn_figure_runner(
            distribution    'norm'      perform a serial transfer
                                        experiment with glenn_epoch_runner(), 
                                        with a set epoch length and transfer 
                                        size.
                            'expd'      perform a serial transfer with
                                        epoch lengths drawn from an
                                        exponential distribution.
                            'tran'      perform a serial transfer with
                                        transfer sizes drawn from a
                                        lognormal distribution.
                            'mevo'      perform a serial transfer where the
                                        regulator and non-regulator lists
                                        are not used to create separate
                                        phenotypes to compete, but instead
                                        are used to create a phenotype
                                        matrix for co-evolution of the
                                        traits.
            eplen           ?           the mean length of an epoch in the
                                        serial transfer experiment
            steps           ?           the size of a single intergation
                                        step.
            maxeps          ?           maximal number of epochs of the
                                        experiment. If the experiment
                                        comes within equilibrium tolerance,
                                        it is cut shorter.
            philist         []          list of lysogeny propensities.
            taulist         []          list of regulation thresholds.
            birth           ?           host birthrate.
            death           ?           host deathrate.
            resurge         ?           prophage induction rate.
            infectivity     ?           virion infectivity.
            burstsize       ?           number of virion progeny per lysis.
            decay           ?           virion decay rate.
            fomut           ?           forward mutation rate.
            bamut           ?           back mutation rate.
            inits           []          list of initial phenotype
                                        concentrations.

            
%}

    tolrnc=1e-7;

    if distribution=='norm'
        [tf,cf,~,~,te,ce,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner(eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,inits,tolrnc);
    else
        if distribution=='expd'
            [tf,cf,~,~,te,ce,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner_exprnd(eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,inits,tolrnc);
        else
            if distribution=='tran'
                [tf,cf,~,~,te,ce,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner_transfer(eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,inits,tolrnc);
            else
                if distribution =='mevo'
                	[tf,cf,~,~,te,ce,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner_morevo(eplen,steps,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,inits,tolrnc);
                end
            end
        end
    end

    ce=cat(2,inits',ce);
%   FOR NON-REGULATORS (NO SIGNAL PLOT)
%   all
    glenn_figure_awesome(tf,cf(1,:),cf(2:1+length(philist)+length(taulist),:),cf(2+length(philist)+length(taulist):end,:),philist,te,ce(1,:),ce(2:1+length(philist)+length(taulist),:),ce(2+length(philist)+length(taulist):end,:));
%   line
    glenn_figure_awesome(tf,cf(1,:),cf(2:1+length(philist)+length(taulist),:),cf(2+length(philist)+length(taulist):end,:),philist,[],[],[],[]);
%   dash
    glenn_figure_awesome([],[],[],[],philist,te,ce(1,:),ce(2:1+length(philist)+length(taulist),:),ce(2+length(philist)+length(taulist):end,:));
         
%   Calculate signal:
sig1=zeros(1,length(tf));
for pp=1:size(cf,2)
    sig1(pp)=.25*sum(cf(2+length(philist)+length(taulist):end,pp))*cf(1,pp);
end
sig2=zeros(1,length(te));
for pp=1:size(ce,2)
    sig2(pp)=.25*sum(ce(2+length(philist)+length(taulist):end,pp))*ce(1,pp);
end

%   FOR REGULATORS (SIGNAL PLOT PRESENT)
%   all
    glenn_figure_awesome_tau(tf,cf(1,:),cf(2:1+length(philist)+length(taulist),:),cf(2+length(philist)+length(taulist):end,:),taulist,te,ce(1,:),ce(2:1+length(philist)+length(taulist),:),ce(2+length(philist)+length(taulist):end,:),sig1,sig2);
%   line
    glenn_figure_awesome_tau(tf,cf(1,:),cf(2:1+length(philist)+length(taulist),:),cf(2+length(philist)+length(taulist):end,:),taulist,[],[],[],[],sig1,[]);
%   dash
    glenn_figure_awesome_tau([],[],[],[],taulist,te,ce(1,:),ce(2:1+length(philist)+length(taulist),:),ce(2+length(philist)+length(taulist):end,:),[],sig2);
          

    
end