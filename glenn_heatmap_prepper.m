function glenn_heatmap_prepper(philist,taulist,folder,phiortau,distribution,identifier,eplist)

% Data treatment and figure making figure to work with glenn_heatmap_runne*() function output.

    %{
        This function can do two things.
            In the case of a heatmap run that performed a co-evolution of
            two parameters, it can turn that into a heatmap of the
            phentype space for each epoch length it is fed.
            In the case of a heatmap run that evolved only one parameter,
            it can produce a single heatmap, containing data from all the
            serial transfer experiments with different epoch lengths or
            transfer sizes.

        Ordinarily, this function is called by other functions, and only
        handles the data processing and graphical presentation!

        Usage:
            glenn_heatmap_prepper(
                    philist         list of maximal lysogeny propensities
                    taulist         list of regulation thresholds
                    folder          folder containing the data files
                    phiortau        first name flag
                    distribution    second name flag, decideds whether tio treat
                                    the input data as co-evolution run or not  
                    identifier      third name flag
                    eplist          list of epochs which are to be included
                                    in the making of the heatmap
                )

        This function reads from files produced by heatmap runner
        functions, and turns them into figures after processing. It uses
        two sub-functions given below.

    %}

    if distribution ~= "mevo"
    	list=[philist taulist];
        num=length(list);
        matrix=[];
        signal=[];
        for eplen = eplist
            run=0;
            fprintf("WORKING ON EPLEN "+num2str(eplen)+'\n')
            fc=fixer(folder+phiortau+"_"+distribution+"_"+identifier+"_"+"run_"+num2str(eplen,'%i')+"_fc"+".txt",num);
            pc=zeros(num,size(fc,2));
            mx=zeros(1,size(fc,2));
            for jj = 1:size(fc,2)
                pc(1:num,jj)=(fc(2:(1+num),jj)+fc((2+num):(1+2*num),jj));
                pc(1:num,jj)=pc(1:num,jj)/sum(pc(1:num,jj));
                mx(1,jj)=0.25*fc(1,jj)*sum(fc((2+num):(1+2*num),jj));
            end
            sigmax=max(mx);
            sigmea=mean(mx);
            sigend=mx(1,jj);
            pca=zeros(num,1);
            for ii = 1:num
                if size(pc,2)>2000*4*eplen
                    pca(ii,1)=mean(pc(ii,(2000*4*eplen):end));
                else
                    run=1;
                    pca(ii,1)=mean(pc(ii,(200*4*eplen):end));
                end
            end
            if run==1
                fprintf("RUN TOO SHORT, TOOK MEAN OVER EPOCH 200 TO END\n")
            end
            if sum(pca(:,1))<1-(1e-5) || sum(pca(:,1))>1+(1e-5)
                fprintf("WRONG SUM IN MATRIX COLUMN: "+num2str(sum(pca(:,1)),'%f')+"\n")
                % 'Cuz numerical errors gone big, yo.
            end
            matrix=cat(2,matrix,pca(:,1));
            signal=cat(2,signal,[sigmax;sigmea;sigend]);
        end
        
        fprintf("WRITING TO FILES"+'\n')
        matrix_file=folder+phiortau+"_"+distribution+"_"+identifier+"_"+"m"+".txt";
        xticklabellist_file=folder+phiortau+"_"+distribution+"_"+identifier+"_"+"x"+".txt";
        yticklabellist_file=folder+phiortau+"_"+distribution+"_"+identifier+"_"+"y"+".txt";
        signal_file=folder+phiortau+"_"+distribution+"_"+identifier+"_"+"s"+".txt";
        
        xaxislabeltext="Distribution of $\"+phiortau+"$ in population";
        yaxislabeltext="Length of epoch in ($b^{-1}$)";
        
        dlmwrite(matrix_file,matrix');
        dlmwrite(xticklabellist_file,list);
        dlmwrite(yticklabellist_file,eplist);
        dlmwrite(signal_file,signal);
        
        glenn_heatmap_maker(matrix_file,'jet',xticklabellist_file,yticklabellist_file,xaxislabeltext,yaxislabeltext,signal_file)

    else
        plist=philist;
        tlist=taulist;
        num=length(philist)*length(taulist);
        for eplen = eplist
            fprintf("WORKING ON EPLEN "+num2str(eplen)+'\n')
            fc=importdata(folder+phiortau+"_"+distribution+"_"+identifier+"_"+"run_"+num2str(eplen,'%i')+"_fc"+".txt");
            pc=zeros(num,size(fc,2));
            for jj = 1:size(fc,2)
                pc(1:num,jj)=(fc(2:(1+num),jj)+fc((2+num):(1+2*num),jj));
                pc(1:num,jj)=pc(1:num,jj)/sum(pc(1:num,jj));
            end
            pca=zeros(num,1);
            for ii = 1:num
                pca(ii,1)=mean(pc(ii,2000:end));
            end
            if sum(pca(:,1))<1-(1e-5) || sum(pca(:,1))>1+(1e-5)
                fprintf("WRONG SUM IN MATRIX COLUMN: "+num2str(sum(pca(:,1)),'%f')+'\n')
                % 'Cuz numerical errors gone big, yo.
            end
            matrix=zeros(length(plist),length(tlist));
            for t = 1:length(tlist)
                for p = 1:length(plist)
                    matrix(t,p)=pca((t-1)*length(plist)+p,1);
                end
            end
            
            fprintf("WRITING TO FILES"+'\n')
            matrix_file=folder+phiortau+"_"+distribution+"_"+identifier+"_"+"m"+"_"+num2str(eplen,'%i')+".txt";
            xticklabellist_file=folder+phiortau+"_"+distribution+"_"+identifier+"_"+"x"+"_"+num2str(eplen,'%i')+".txt";
            yticklabellist_file=folder+phiortau+"_"+distribution+"_"+identifier+"_"+"y"+"_"+num2str(eplen,'%i')+".txt";
            
            xaxislabeltext="Distribution of $\phi_{max}$ in population";
            yaxislabeltext="Distribution of $\tau$ in population";
            
            dlmwrite(matrix_file,matrix);
            dlmwrite(xticklabellist_file,plist);
            dlmwrite(yticklabellist_file,tlist);
            
            glenn_heatmap_maker(matrix_file,'pink',xticklabellist_file,yticklabellist_file,xaxislabeltext,yaxislabeltext,[])
            
        end
        
    end
    
    fprintf("DONE"+'\n')
    
end

% UGLY SUB-FUNCTION TO PREVENT SHIT FROM GOING WEIRD
% A rare occurrence where the file written before is not properly newlined,
% this function ought to fix that and make it properly readable for this
% file.
function matfixed=fixer(naam,list)
    types=list*2+1;
    fprintf("READING FILE...\n")
    mat=importdata(naam);
    if size(mat,1)==1
        fprintf(" FILE BROKEN, FIXING...")
        matfixed=zeros(types,size(mat,2)/types);
        for ii = 1:types
            matfixed(ii,:)=mat(((ii-1)*size(mat,2)/types+1):(ii*size(mat,2)/types));
        end
        fprintf(" -> FIXED, CONTINUE\n")
    else
        if size(mat,1)~=types
            fprintf("\nDAFUQ?!\n")
        else
            matfixed=mat;
            fprintf(" ALL FINE, CONTINUE\n")
        end
    end
end


function glenn_heatmap_maker(matrix_file,~,xticklabellist_file,yticklabellist_file,yaxislabeltext,xaxislabeltext,signal_file)
    
    matrix=importdata(matrix_file)';
    yticklabellist=importdata(xticklabellist_file)';
    xticklabellist=importdata(yticklabellist_file);
    if isempty(signal_file)==0
        signal=importdata(signal_file)';
    else
        signal=[];
    end
    
    figure;
    host=axes('position',[.1 .1 .85 .85]);
    set(host,'TickLabelInterpreter','latex');
    
    colormap(flipud(bone));
    map=imagesc(matrix);
    set(map,'Parent', host)
    set(host,'Ydir','normal')
    
    xticklist=1:length(xticklabellist);
    yticklist=1:length(yticklabellist);
    
    xticks(xticklist)
    xticklabels(xticklabellist)
    ylabel(xaxislabeltext)
    

    if isempty(signal)==0
        hold on
        
        conversion=yticklist(end)/yticklabellist(end);
        
        line1=plot(xticklist,signal(:,1).*conversion(end),'Color','r','Marker','none','LineWidth',1,'linestyle',':');
        line2=plot(xticklist,signal(:,2).*conversion(end),'Color','r','Marker','none','LineWidth',1,'linestyle',':');
        line3=plot(xticklist,signal(:,3).*conversion(end),'Color','r','Marker','none','LineWidth',2);
        set(line1,'Parent',host)
        set(line2,'Parent',host)
        set(line3,'Parent',host)
        
    end    

    yticks(yticklist)
    yticklabels(yticklabellist)
    xlabel(yaxislabeltext)
    colorbar;
    set(host,'TickLabelInterpreter','latex');

    
end
