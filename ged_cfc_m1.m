function [mvarpac,zTOT,evecsT,evals,epochs,comps,hicomps,ERRORS]=ged_cfc_m1(X,chanlocs,sr,lofreq,hifreq,lo_locked_TF_plots,skip_pac,comps,hicomps,zTOT)

%% [mvarpac,zTOT,evecsT,evals,epochs,comps,hicomps,ERRORS]=ged_cfc_m1(X,chanlocs,sr,lofreq,hifreq,lo_locked_TF_plots,savefn,comps,hicomps,zTOT)
% This function computes multivariate cross-frequency coupling, based on 
% generalized eigendecomposition (gedCFC; Cohen, 2017). This uses method
% one from the Cohen paper and focuses on trough-locked activity. It is
% currently minimally flexible and requires most inputs. 
%
%  Inputs:
%     X -- EEG data, channels X timepoints, any trial data will be vectorized
%     chanlocs -- channel locations structure from EEGlab
%     sr -- sampling rate [5000]
%     lofreq -- low frequency for component selection [peak alpha frequency]
%     hifreq -- high frequency range for PAC computation [18:80]
%     lo_locked_TF_plots -- if true, trough-locked TF plots are calculated and plotted based on components selected. [false]
%     savefn -- not functional yet, does nothing
%     L -- pca loadings if pca activations is passed
%   
%   Time-saving Inputs/Outputs: In the event of errors, these will be returned if
%   requested as outputs. These outputs can be used as incomplete inputs to
%   the function to speed up the process and prevent unnnecessary
%   repetition of steps.
%     comps -- vector of low frequency components
%     hicomps -- cell array {1 X numel(comps)} of vectors for high frequency components
%     zTOT -- cell array {1 X numel(comps)} of trough-locked TF plots for use with making Canolty plots.
% 
%  Outputs:
%     mvarpac -- multivariate phase-amp coupling estimates. [hicomps X freq X locomps], data is zero where the number of hi comps differs across low comps
%     zTOT -- see above,,,, 
%     evecsT -- 
%     evals -- 
%     epochs -- cell array {1 X numel(comps)} for trough-locked indices
%     comps -- see above
%     hicomps -- see above
%     ERRORS -- call rethrow(ERRORS) to return any errors that prevented the function from executing.

% Method 1: covariance matrix at a low-frequency phase value, 
% relative to covariance from entire time period.

% You will need the following files in the current directory or Matlab path:
%   - emptyEEG.mat 
%   - filterFGx.m
%   - topoplotIndie.m

% mikexcohen@gmail.com

%% initialize variables
% if ~exist('use_pca','var') || isempty(use_pca)
%     use_pca=false;
% else
%     [L,PC,E]=princomp(X(:,nnans)'); 
%     E=E/sum(E); 
%     varexp=round(sum(E(1:2))*1000)/10;
%     ncomp=find(varexp>0,1,'last');
%     L=L(:,1:ncomp); PC=PC(:,1:ncomp)'; E=E(1:ncomp);
%     X=reshape(PC(:,1:ncomp)',[ncomp size(X,2) size(X,3)]);
% end

if nargin==0; 
    display('[mvarpac,zTOT,evecsT,evals,epochs,comps,hicomps,ERRORS]=ged_cfc_m1(X,chanlocs,sr,lofreq,hifreq,lo_locked_TF_plots,skip_pac,comps,hicomps,zTOT)')
    return
end

if ~exist('lofreq','var') || isempty(lofreq)
    display('No low frequency range given, searching for peak frequency.')
    ch1=~any(isnan(X(:,:)),2); bch1=find(~ch1);
    [PSD,f]=pwelch(X(ch1,:)',sr,0,[],sr,'psd');
    PSD=mean(bsxfun(@rdivide,PSD,sum(PSD,1)),2);
    PSD=PSD(f>0&f<50,1); f=f(f>0&f<50);
        [TREND]=polyfit(log10(f),log10(PSD),3);
        Y=polyval(TREND,log10(f));
        lopeak=find((log10(PSD)-Y)==max(log10(PSD)-Y));
            figure; plot(log10(f),log10(PSD));
            hold on;
            plot(log10(f),Y);
            plot(log10(f),log10(PSD)-Y);
            plot(log10(f(lopeak)),log10(PSD(lopeak))-Y(lopeak),'o');
        alpha=f(lopeak);
        lofreq=round(alpha,1);
        display(['    Continuing with alpha peak identified at ' num2str(round(alpha,1)) '.'])
end
if ~exist('hifreq','var') || isempty(hifreq)
    hifreq = 18:80;
end
if ~exist('lo_locked_TF_plots','var') || isempty(hifreq)
    lo_locked_TF_plots=false;
end
TIME=tic;

try
%% GED to identify theta component
% theta covariance
lofilt = filterFGx(X(:,:),sr,lofreq,5);
lofilt = bsxfun(@minus,lofilt,mean(lofilt,2));
locov  = (lofilt*lofilt')/size(X,2);

% broadband covariance
tmpdat = bsxfun(@minus,X(:,:),mean(X(:,:),2));
bbcov  = (tmpdat*tmpdat')/size(X,2);

% GED
[evecsT,evals] = eig(locov,bbcov);

% find best component and compute filter projection
[~,maxcomp] = sort(diag(evals));
lomap    = locov*evecsT;
%%% for future pca variant, use pinv(L(:,a:b))*locov*evecsT;
% lomap    = lomap(:,maxcomp(end));
% lomap    = lomap(:,maxcomp);

% fix sign of map (max is positive)
[~,maxe] = max(abs(lomap));
lomap = lomap .* repmat(sign(lomap(maxe)),[size(lomap,1) 1]);

% theta time series component
locomp = lofilt' * evecsT(:,1); % do the first one
locomp = locomp * sign(corr(locomp,filterFGx(X(maxe(1),:),sr,lofreq,9)'));
 for i=2:size(evecsT,2);
     locomp(:,i) = lofilt' * evecsT(:,i); % now loop through the rest of them
     % fix sign of time series according to sign of correlation with EEG
     locomp(:,i) = locomp(:,i) * sign(corr(locomp(:,i),filterFGx(X(maxe(i),:),sr,lofreq,9)'));
 end
 display(['Low frequency components isolated in ', num2str(round(toc(TIME(1)))), ' secs.'])

%% identify troughs and get surrounding covariance matrices
nwin = ceil(sr/lofreq/8); % window size is 1/4 cycle (1/8 of either side)
if ~exist('comps','var') || isempty(comps)
    figure(10); clf
    for i=1:size(locomp,2)
        subplot(5,7,i);
        topoplotIndie(abs(lomap(:,i)), chanlocs,'numcontour',0,'electrodes','on','shading','interp'); 
        title(['Comp ' int2str(i)]);
    end
    prompt = 'Select low frequency components of interest...';
    comps=[input(prompt)];
        if isempty(comps);
            display('  ...  continuing with all components')
            comps=1:size(locomp,2);
        else
            display(['  ugh! finnnnee, continuing with components ' int2str(comps)])
        end
end
% find troughs
troughs={}; epochs={};
for i=1:length(comps); 
    [~,trawfs]=findpeaks(-locomp(:,comps(i)),'MinPeakDistance',nwin);
    trawfs(trawfs<sr) = [];    %troughs(troughs<nwin+1) = [];
    trawfs(trawfs>size(X(:,:),2)-sr) = [];    %troughs(troughs>size(X(:,:),2)-nwin-1) = [];
        troughs{i}=trawfs;
        segs=[ones(length(trawfs),1) trawfs-2500 trawfs+2500];
        segs=remove_edge_epochs(segs,12500); % this line is hard coded for impa1 resting state data (1.25 sec segments)
        epochs{i}=segs;
end
        
covT = zeros(size(X,1),size(X,1),length(troughs));
evecs=zeros(size(covT)); evals=evecs; himaps=evecs;

% trough-locked covariance
for ci=1:length(comps)
    for ti=1:length(troughs{ci})
        tmpdat = X(:,(troughs{ci}(ti)-nwin):(troughs{ci}(ti)+nwin));
        tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
        covT(:,:,ci) = covT(:,:,ci) + (tmpdat*tmpdat')/nwin;
    end
    covT(:,:,ci) = covT(:,:,ci)./ti;

    %%% GED to get gamma peak/trough networks %%%
    [evecs(:,:,ci),evals(:,:,ci)] = eig(covT(:,:,ci),bbcov);
    [~,compidx]   = sort(diag(evals(:,:,ci))); % max component
        evecs(:,:,ci)=evecs(:,compidx,ci);
        evals(:,:,ci)=evals(:,compidx,ci);
        
    himaps(:,:,ci) = covT(:,:,ci)*evecs(:,:,ci); % forward model of filter
end

if ~exist('hicomps','var') || isempty(hicomps)
    hicomps={}; ask_comps=true;
    else ask_comps=false;
end
for ci=1:length(troughs)

    if ask_comps
    figure(10); clf 
    for i=1:size(locomp,2)
        subplot(5,7,i);
        topoplotIndie(abs(himaps(:,i,ci)), chanlocs,'numcontour',0,'electrodes','on','shading','interp'); 
        title(['Comp ' int2str(i)]);
    end

    prompt = ['Select high frequency components of interest for component ' int2str(comps(ci)), '...'];
    select=[input(prompt)];
        if isempty(select);
            display('  ...  continuing with all components')
            select=1:size(locomp,2);
            hicomps{ci}=select;
        else
            hicomps{ci}=select;
            display(['  I guess that''s cool. Continuing with components ' int2str(hicomps{ci})])
        end
    end



    hinet = himaps(:,hicomps{ci},ci);
    % fix sign
    [~,idx] = max(abs(hinet));
    for i=1:length(idx); 
        hinet(:,i) = hinet(:,i) * sign(hinet(idx(i),i));
    end
end

maxcomps= max(cellfun('length',hicomps));
    
%% get time course and reconstruct topography and sources
if ~skip_pac

display('Beginning PAC calculations...')
for ci=1:length(comps)
for hci=1:length(hicomps{ci})
    
    if ci==1 && hci==1; %% initialize mvarpac on first 
        mvarpac = zeros(maxcomps,length(hifreq),length(comps));
    end

    parfor fi=1:length(hifreq)

        % bandpass filter trough component
        hicomp = filterFGx(X(:,:),sr,hifreq(fi),15)' * evecs(:,hicomps{ci}(hci),ci); %#ok<PFTUS>
        hicomp=hicomp./max(abs(hicomp(:)));
        % find peaks and troughs and get distances
        hiwin = ceil(sr/hifreq(fi)/8);
        [~,troughs] = findpeaks(gaussmooth_signal(-hicomp,10),'MinPeakDistance',hiwin);
        [~,peeks] = findpeaks(gaussmooth_signal(hicomp,10),'MinPeakDistance',hiwin);   
        mvarpac(hci,fi,ci) = mean(hicomp(peeks)) - mean(hicomp(troughs));
    end
end
    display(['PAC computed for one low frequency component at ', num2str(round(toc(TIME(1))/60,2)), ' min.'])

    if lo_locked_TF_plots

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



end

ERRORS=[];

display(['Total process completed in ', num2str(round(toc(TIME(1))/60,2)), ' minutes'])
display(['    Avg Time Per Component = ', num2str(round(toc(TIME(1))/60,2)/length(comps)), ' minutes'])
% 
% if exist('savefn','var') || ~isempty(savefn)
%     
% end
end
catch ME
    if ~exist('mvarpac','var'); mvarpac=[]; end
    if ~exist('zTOT','var'); zTOT=[]; end
    if ~exist('evecsT','var'); evecsT=[]; end
    if ~exist('epochs','var'); epochs=[]; end
    if ~exist('comps','var'); comps=[]; end
    if ~exist('hicomps','var'); hicomps=[]; end
%     keyboard
    ERRORS=ME;

end


