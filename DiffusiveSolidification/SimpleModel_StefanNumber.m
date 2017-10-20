%% Simple Model Stefan Number Investigation
%
% This script follows on from the mushy factor investigation by also
% considering salt diffusion in the mushy phase. We now use the Stefan
% number rather than the mushy diffusivity factor as our independent
% variable since the mushy diffusivity factor varies with the Lewis number.
%
% It creates figure 3.4 of my pre-corrections thesis.
%
% (02/03/16)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Settings

% Investigation settings
StefanVec   = [0,10.^(-2:0.1:3)];
LqdusInfVec = 0.01:0.01:0.99;
LewisVec    = 10.^(0.6:0.2:5);
Threshold   = 0.01;


%% Main code

% Initialise results container
GrowthRate     = NaN(length(StefanVec),length(LqdusInfVec));
ThresholdLewis = NaN(length(StefanVec),length(LqdusInfVec));
GrowthRate_WS_100 = NaN(length(StefanVec),length(LqdusInfVec));

% Loop over settings
tic
for iLi = 1:length(LqdusInfVec)
    parfor iSt = 1:length(StefanVec)
        
        % Initialise approximate threshold
        ApprxThreshold = NaN;
        ApprxGrowthRate = NaN;
        
        % Calculate growth rate (no salt)
        MushDiff = sqrt(1+StefanVec(iSt));
        PAR_GrowthRate = DiffGrow_Simp_Cond(MushDiff,LqdusInfVec(iLi)); %#ok<PFBNS>
        GrowthRate(iSt,iLi) = PAR_GrowthRate;
        
        % Loop over salt diffusion in liquid phase
        Searching = true;
        GrowthRate_tmp = NaN(size(LewisVec));
        for iS = length(LewisVec):-1:1
            
            % Set initial guess
            if iS == length(LewisVec)
                IG = PAR_GrowthRate;
            else
                IG = GrowthRate_tmp(iS+1);
            end
            
            % Calculate growth rate (with salt)
            MushDiff = sqrt((1+StefanVec(iSt))/(1+StefanVec(iSt)/LewisVec(iS)));
            GrowthRate_tmp(iS) = DiffGrow_Simp_Cond(MushDiff,LqdusInfVec(iLi),sqrt(LewisVec(iS)),IG);
            
            % Check numerical results
            if isnan(GrowthRate_tmp(iS)) || GrowthRate_tmp(iS) < 0
                GrowthRate_tmp(iS) = NaN; %#ok<NASGU>
                break
            end
            
            % Get threshold Lewis number
            if (RelErr(GrowthRate_tmp(iS),PAR_GrowthRate,1) < Threshold) && Searching
                ApprxThreshold = LewisVec(iS);
                ApprxGrowthRate = GrowthRate_tmp(iS);
            else
                Searching = false;
            end
        end
        
        % Salt effects
        try
            
            % Threshold Lewis number
            MushDiff = @(Lewis) sqrt((1+StefanVec(iSt))/(1+StefanVec(iSt)/Lewis));
            GrowthRate_WS_Fcn = @(Lewis) DiffGrow_Simp_Cond(MushDiff(Lewis),LqdusInfVec(iLi),sqrt(Lewis),ApprxGrowthRate);
            ErrFcn = @(Lewis) RelErr(GrowthRate_WS_Fcn(Lewis),PAR_GrowthRate,1)-Threshold;
            ThresholdLewis(iSt,iLi) = fzero(ErrFcn,ApprxThreshold);
            
            % Growth rate for Le = 100
            GrowthRate_WS_100(iSt,iLi) = GrowthRate_WS_Fcn(1000);
        catch
            ThresholdLewis(iSt,iLi) = NaN;
        end
    end
    toc
end

% Tidy
clear iLi iSt ApprxThreshold MushDiff PAR_GrowthRate GrowthRate_tmp IG


%% Calculate effect sign change Stefan number

% Loop over initial liquidus temperatures
SignChangeStefan = NaN(size(LqdusInfVec));
for iLi = 1:length(LqdusInfVec)
    
    % Interpolate data
    LogicalIndex = ~isnan(GrowthRate_WS_100(:,iLi));
    if sum(LogicalIndex) > 3
        SignChangeStefan(iLi) = interp1(GrowthRate(LogicalIndex,iLi)-GrowthRate_WS_100(LogicalIndex,iLi),StefanVec,0);
    end
end

% Tidy
clear iLi LogicalIndex


%% Output
FC = 0;
LW = 3;
FS = 18;

% Conducting growth rates
if 0

    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Graph options
    CTick = [0.1,1,10];
    
    % Create graph
    contourf(StefanVec(2:end),LqdusInfVec,log10(GrowthRate(2:end,:))',50,'edgecolor','none')
    patch(StefanVec([2,2,end,end]),[0,LqdusInfVec( 1 ),LqdusInfVec( 1 ),0],'k')
    patch(StefanVec([2,2,end,end]),[1,LqdusInfVec(end),LqdusInfVec(end),1],'k')
    xlabel('Compositional Stefan number,  S_c','fontsize',FS)
    set(gca,'xscale','log','xtick',10.^(-2:3))
    ylabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS)
    set(gca,'ylim',[0,1],'ytick',0:0.1:1)
    set(gca,'fontsize',FS)
    caxis(log10(CTick([1,end])))
    colormap(cbrewer('seq','YlOrRd',50))
    hC = colorbar('ytick',log10(CTick),'yticklabel',CTick,'fontsize',FS);
    ylabel(hC,'Growth rate,  \lambda','fontsize',FS)
    box on
    
    % Tidy
    clear CTick hC
end

% Threshold Lewis values - Liquid phase only
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Graph options
    CTick = [30,50,100,190,500,1000,2000,5000,10000,20000,50000];
    
    % Create graph
    hold on
    contourf(StefanVec(2:end),LqdusInfVec,log10(ThresholdLewis(2:end,:))',50,'edgecolor','none');
    contour( StefanVec(2:end),LqdusInfVec,log10(ThresholdLewis(2:end,:))',log10([1000,170,100]),'edgecolor','k','linewidth',LW)
    patch(StefanVec([2,2,end,end]),[0,0.03,0.03,0],0.5*[1,1,1])
    patch(StefanVec([2,2,end,end]),[1,LqdusInfVec(end),LqdusInfVec(end),1],0.5*[1,1,1])
    plot(SignChangeStefan,LqdusInfVec,'-','color',[0,0.5,0],'linewidth',LW)
    hold off
    if 0
        xlabel('Compositional Stefan number,  S_c','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(-2:3))
        ylabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS)
        set(gca,'ylim',[0,1],'ytick',0:0.1:1)
        caxis(log10(CTick([1,end])))
        colormap(cbrewer('seq','YlOrRd',50))
        hC = colorbar('ytick',log10(CTick),'yticklabel',CTick,'fontsize',FS);
        ylabel(hC,'Threshold Lewis number,  L_e','fontsize',FS)
        line('parent',hC,'xdata',[-5,5],'ydata',log10( 100*[1,1]),'color','k','linewidth',LW)
        line('parent',hC,'xdata',[-5,5],'ydata',log10( 190*[1,1]),'color','k','linewidth',LW)
        line('parent',hC,'xdata',[-5,5],'ydata',log10(1000*[1,1]),'color','k','linewidth',LW)
        set(gca,'fontsize',FS)
        box on
    else
        xlabel('XLAB','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(-2:3),'xticklabel',{'X01','X02','X03','X04','X05','X06'})
        ylabel('YLAB','fontsize',FS)
        set(gca,'ylim',[0,1],'ytick',0:0.1:1,'yticklabel',{'Y01','Y02','Y03','Y04','Y05','Y06','Y07','Y08','Y09','Y10','Y11'})
        caxis(log10(CTick([1,end])))
        colormap(cbrewer('seq','YlOrRd',50))
        hC = colorbar('ytick',log10(CTick),'fontsize',FS,'yticklabel',{'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12'});
        ylabel(hC,'CLAB','fontsize',FS)
        line('parent',hC,'xdata',[-5,5],'ydata',log10( 100*[1,1]),'color','k','linewidth',LW)
        line('parent',hC,'xdata',[-5,5],'ydata',log10( 190*[1,1]),'color','k','linewidth',LW)
        line('parent',hC,'xdata',[-5,5],'ydata',log10(1000*[1,1]),'color','k','linewidth',LW)
        set(gca,'fontsize',FS)
        box on
    end
    axis square
    
    % Tidy
    clear hC CTick
end

% Tidy
clear FC FS LW

