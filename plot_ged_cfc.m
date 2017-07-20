function [h,zTOT]=plot_ged_cfc(X,hicomps,locomp,evecs,epochs,zTOT)

for ci=1:length(hicomps)
    % recompute hicomps to use for plotting
    for i=1:length(hicomps{ci});
        if i==1;
            hicomp = X(:,:)' * evecs(:,hicomps{ci}(i));
            if ~isreal(hicomp(1,1)); hicomp=abs(hicomp)-mean(abs(hicomp(:,1)),1); end
        else
            hicomp(:,i) = X(:,:)' * evecs(:,hicomps{ci}(i));
            if ~isreal(hicomp(1,i)); hicomp(:,i)=abs(hicomp(:,1))-mean(abs(hicomp(:,1)),1); end
        end
    end
    [v_trough]=epochs2vect(epochs{ci});
    LOW=reshape(locomp(v_trough,comps(ci)),[1 epochs{ci}(1,3)-epochs{ci}(1,2)+1 length(epochs{ci}(:,1))]);
    if ~exist('zTOT','var') || isempty(zTOT{ci})
        TIME(2)=tic;
        HIGH=reshape(hicomp(v_trough,:)',[size(hicomp,2) epochs{ci}(1,3)-epochs{ci}(1,2)+1 length(epochs{ci}(:,1))]);
        [in,ev]=timefq(HIGH,hifreq,sr);
        tot=ev+in; clear in ev
        totz=bsxfun(@rdivide,bsxfun(@minus,tot,nanmean(tot,2)),nanstd(tot,[],2)); clear tot
            zTOT{ci}=totz;
            display(['   TF transform completed in ' num2str(round(toc(TIME(2))/60,2)), ' minutes.'])
    end
    if ci==1; zTOT{numel(comps)}=[]; end    

    figure('units','normalized','outerposition',[0 0 1 1]); c=0; d=1;
        for i=1:length(hicomps{ci});
            if i>1 && mod(i,6)==1;
                figure('units','normalized','outerposition',[0 0 1 1]); c=-6; d=1;
            end
            subplot(4,3,i+c); a=get(gca,'Position');  % left bottom width height
            subplot(4,3,i+3+c); b=get(gca,'Position');  % left bottom width height
            if i==1 || i==4;
                a(1) = a(1) - .06;
                b(1) = b(1) - .06;
            else if i==3 || i==6;
                    a(1) = a(1) + .06;
                    b(1) = b(1) + .06;
                end
            end

            subplot('Position',[a(1)+.055 a(2) a(3)-.055 a(4)]); 
                imagesc(-500:.2:500,hifreq,zTOT{ci}(:,:,i),[-3 3]); 
                axis xy; set(gca,'XTickLabel',[]);
            subplot('Position',[b(1)+.055 b(2) b(3)-.055 b(4)]); 
                plot(-500:.2:500,mean(LOW,3));

            axes('Position',[a(1)-.09 a(2)-(.5*(d-a(2))) .17 .17]); 
                topoplotIndie(abs(himaps(:,hicomps{ci}(i),ci)),chanlocs,'numcontour',0,'electrodes','off','shading','interp');
            axes('Position',[a(1)+a(3) a(2) .03 a(4)]); 
                plot(mvarpac(i,:,ci),hifreq); xlim([min(mvarpac(i,:,ci)) max(mvarpac(i,:,ci))]); ylim([min(hifreq) max(hifreq)]);
                set(gca,'Visible','off');
            if i==3; c=c+3; d=.5; end
            drawnow
        end


end