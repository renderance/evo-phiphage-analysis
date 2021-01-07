function glenn_signal_surface_runner(which,filename,timename,startsname)
% Surface plot generator for signal of regulators and non-regulators using glenn_epoch_runner().
    %{
     glenn_signal_surface_runner() 
     This function performs an invasion experiment of a particular
     phenotype onto a senstive bacterial population, then calculates the
     signal that it results in over time. It does so for multiple
     phenotypes and plots the result as a surface.

       Usage: glenn_solver(
               which = "tau" or "phi", determines whether a signal surface
               plot for regulator or non-regulator phenotypes will be
               constructed.
               filename = name of output file for the signal matrix
               timename = name of output file for the progression of time
               startsname = name of output file with phenotype values to
               start with.
           )

     Note that this function calls the epoch runner script to do its work.

    %}

    % Set all necessary run parameters:
    stepsize = 0.25;
    maxeps = 1;

    % Set the bacteriophage/host parameters:
    birth=1;
    death=.3;
    resurge=.4;
    infectivity=.25;
    burstsize=4;
    decay=.5;
    fomut=1e-3;
    bamut=1e-3;

    % Determine the time for which signal should be plotted:
    eplen=125;
    
    % Create an empty signal object.
    signal=[];
    
    if which == "tau"
        % Create a list of regulator phenotypes.
        starts=linspace(0,0.11,100);
        philist=[];
    else
        if which == "phi"
            % Create a list of non-regulator phenotypes.
            taulist=[];
            starts=linspace(0,1,100);
        end
    end
    
    % Go over all the start values, and put them in the right list.
    for ii = starts
        if which=="tau"
            taulist=ii;
        else
            if which=="phi"
                philist=ii;
            end
        end
        
        % Then create a starting conditions vector, which should hold 3 values.
        ins(1)=.7;
        ins(2)=0;
        ins(3)=0.001;
        
        % Use the epoch runner script to solve the invasion process for
        % the bacteriophage phenotype under study.
        [tf,cf,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = glenn_epoch_runner(eplen,stepsize,maxeps,philist,taulist,birth,death,resurge,infectivity,burstsize,decay,fomut,bamut,ins,0);
        
        % Make empty vector for the total bacteriophage population.
        phagesum=zeros(1,size(cf,2));
        % And for the sensitive population.
        sens=zeros(1,size(cf,2));
        % And for the signal.
        sig=zeros(1,size(cf,2));
        
        % Go over the full outcome matrix and calculate sensitives, phages
        % and signal values. Put those in the empty vectors.
        for jj = 1:size(cf,2)
            phagesum(jj)=sum(cf(1+length(philist)+length(taulist)+1:end,jj));
            sens(jj)=cf(1,jj);
            sig(jj)=infectivity*phagesum(jj)*sens(jj);
        end
        
        % Append the outcome of this run to the signal value.
        signal=cat(1,signal,sig);
    end
    
    % Write the output to files.
    dlmwrite(filename,signal);
    dlmwrite(timename,tf);
    dlmwrite(startsname,starts);
    
    % And read those files to create a beautiful image.
    glenn_signal_surface_maker(timename,startsname,filename,"Time","Phenotypes");

    
    function glenn_signal_surface_maker(x_file,y_file,m_file,x_text,y_text)
        x=importdata(x_file);
        y=importdata(y_file);
        m=importdata(m_file);
        figure;
        surf(x,y,m,'LineStyle',':');
        xlabel(x_text);
        xlim([0 x(end)]);
        ylabel(y_text);
        ylim([0 y(end)]);
        zlabel('Signal');
        zlim([0 0.11]);
    end
    
    
end


