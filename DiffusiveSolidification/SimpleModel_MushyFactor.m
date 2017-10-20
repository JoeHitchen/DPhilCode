%% Simple Model Mushy Factor Investigation
%
% This script investigates the effect of the mushy diffusivity factor and
% initial liquidus temperature on ice growth rates in the simplifed model.
% It also investigates the Lewis number needed in the liquid phase to cause
% a substantial change in the growth rate.
%
% It creates figures 3.1, 3.2, and 3.3 of my pre-corrections thesis.
%
% (02/03/16)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Settings

% Investigation vectors
MushyDiffVec = 10.^(0:0.025:3);
LqdusInfVec  = 0.01:0.01:0.99;

% Calculate Salt Effects
SaltEffects = true;
LewisVec    = 10.^(1:0.5:6);
Threshold   = 0.01;


%% Main code

% Initialise results container
GrowthRate      = NaN(length(MushyDiffVec),length(LqdusInfVec));
GrowthRate_Aprx = NaN(length(MushyDiffVec),length(LqdusInfVec));
GrowthRate_Emms = NaN(length(MushyDiffVec),length(LqdusInfVec));
if SaltEffects
    ThresholdLewis = NaN(length(MushyDiffVec),length(LqdusInfVec));
end

% Loop over settings
tic
for iLi = 1:length(LqdusInfVec)
    parfor iM = 1:length(MushyDiffVec)
        
        % Initialise approximate threshold
        ApprxThreshold = NaN;
        
        % Calculate growth rate (no salt)
        PAR_GrowthRate = DiffGrow_Simp_Cond(MushyDiffVec(iM),LqdusInfVec(iLi)); %#ok<PFBNS>
        GrowthRate(iM,iLi) = PAR_GrowthRate;
        
        % Approximate growth rates
        GrowthRate_Aprx(iM,iLi) = 2/MushyDiffVec(iM)*sqrt(log(1/(1/LqdusInfVec(iLi)-1)*MushyDiffVec(iM)));
        GrowthRate_Emms(iM,iLi) = 2/MushyDiffVec(iM)*sqrt(log(1/(1/LqdusInfVec(iLi)-1)));
        if ~isreal(GrowthRate_Aprx(iM,iLi)); GrowthRate_Aprx(iM,iLi) = NaN; end
        if ~isreal(GrowthRate_Emms(iM,iLi)); GrowthRate_Emms(iM,iLi) = NaN; end
        
        % Calculate salt effects
        if SaltEffects
    
            % Loop over salt diffusion in liquid phase
            GrowthRate_tmp = NaN(size(LewisVec));
            for iS = length(LewisVec):-1:1
                
                % Set initial guess
                if iS == length(LewisVec)
                    IG = PAR_GrowthRate;
                else
                    IG = GrowthRate_tmp(iS+1);
                end
                
                % Calculate growth rate (with salt)
                GrowthRate_tmp(iS) = DiffGrow_Simp_Cond(MushyDiffVec(iM),LqdusInfVec(iLi),sqrt(LewisVec(iS)),IG);
        
                % Check numerical results
                if isnan(GrowthRate_tmp(iS)) || GrowthRate_tmp(iS) < 0
                    GrowthRate_tmp(iS) = NaN; %#ok<NASGU>
                    break
                end
            
                % Get threshold Lewis number
                if RelErr(GrowthRate_tmp(iS),PAR_GrowthRate,1) < Threshold
                    ApprxThreshold = LewisVec(iS);
                end
            end
    
            % Salt effects
            try
                GrowthRate_WS_Fcn = @(Lewis) DiffGrow_Simp_Cond(MushyDiffVec(iM),LqdusInfVec(iLi),sqrt(Lewis),PAR_GrowthRate);
                ErrFcn = @(Lewis) RelErr(GrowthRate_WS_Fcn(Lewis),PAR_GrowthRate,1)-Threshold;
                ThresholdLewis(iM,iLi) = fzero(ErrFcn,ApprxThreshold);
            catch
                ThresholdLewis(iM,iLi) = NaN;
            end
        end
    end
    toc
end

% Tidy
clear iLi iM ApprxThreshold PAR_GrowthRate iS IG GrowthRate_tmp GrowthRate_WS_Fcn ErrFcn


%% Output
FC = 0;
LW = 3;
FS = 18;

% Conducting growth rates
if 1

    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Graph options
    %CTick = [0.004,0.007,0.01,0.02,0.04,0.07,0.1,0.2,0.4,0.7,1,2];
    CTick = [0.003,0.01,0.032,0.100,0.32,1,3.2];
    
    % Create graph
    hold all
    contourf(MushyDiffVec,LqdusInfVec,log10(GrowthRate)',50,'edgecolor','none')
    contour(MushyDiffVec,LqdusInfVec,log10(GrowthRate)',log10(CTick),'edgecolor','k','linewidth',LW)
    patch(MushyDiffVec([1,1,end,end]),[0,LqdusInfVec( 1 ),LqdusInfVec( 1 ),0],0.5*[1,1,1])
    patch(MushyDiffVec([1,1,end,end]),[1,LqdusInfVec(end),LqdusInfVec(end),1],0.5*[1,1,1])
    hold off
    if 0
        xlabel('Mushy diffusivity factor,  m','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(0:3))
        ylabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS)
        set(gca,'ylim',[0,1],'ytick',0:0.1:1)
        set(gca,'fontsize',FS)
        colormap(cbrewer('seq','YlOrRd',50))
        caxis(log10([min(CTick),max(CTick)]))
        hC = colorbar('ytick',log10(CTick),'yticklabel',CTick,'fontsize',FS);
        ylabel(hC,'Growth rate,  \lambda','fontsize',FS)
        for iC = CTick
            line('parent',hC,'xdata',[-5,5],'ydata',log10(iC)*[1,1],'linewidth',LW)
        end
        box on
    else
        xlabel('XLAB','fontsize',FS)
        set(gca, 'xscale', 'log', 'xtick', 10.^(0:3), 'xticklabel', {'X1','X2','X3','X4'})
        ylabel('YLAB','fontsize',FS)
        set(gca,'ylim',[0,1],'ytick',0:0.1:1, 'yticklabel', {'Y01','Y02','Y03','Y04','Y05','Y06','Y07','Y08','Y09','Y10','Y11'})
        set(gca,'fontsize',FS)
        colormap(cbrewer('seq','YlOrRd',50))
        caxis(log10([min(CTick),max(CTick)]))
        hC = colorbar('ytick',log10(CTick),'yticklabel',{'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12'},'fontsize',FS);
        for iC = CTick
            line('parent',hC,'xdata',[-5,5],'ydata',log10(iC)*[1,1],'linewidth',LW)
        end
        ylabel(hC,'CLAB','fontsize',FS)
        box on
    end
    axis square
    
    % Tidy
    clear CTick hC
end


% Approximate growth rates
if 1
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Graph options
    CTick = [0.001,0.002,0.004,0.006,0.01,0.02,0.04,0.06,0.1,0.2,0.4,0.6,1];
    
    % Create graph
    hold on
    contourf(MushyDiffVec,LqdusInfVec,log10(RelErr(GrowthRate_Aprx,GrowthRate,1))',50,'edgecolor','none')
    contour( MushyDiffVec,LqdusInfVec,log10(RelErr(GrowthRate_Aprx,GrowthRate,1))',log10([0.001,0.01,0.1]),'edgecolor','k','linewidth',LW)
    plot(MushyDiffVec,1./(MushyDiffVec+1),'color',[0,0.5,0],'linewidth',LW)
    hold off
    patch(MushyDiffVec([1,1,end,end]),[0,LqdusInfVec( 1 ),LqdusInfVec( 1 ),0],0.5*[1,1,1])
    patch(MushyDiffVec([1,1,end,end]),[1,LqdusInfVec(end),LqdusInfVec(end),1],0.5*[1,1,1])
    if 0
        xlabel('Mushy diffusivity factor,  m','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(0:3))
        ylabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS)
        set(gca,'ylim',[0,1],'ytick',0:0.1:1)
        caxis(log10(CTick([1,end])))
        colormap(cbrewer('seq','YlOrRd',50))
        hC = colorbar('ytick',log10(CTick),'yticklabel',100*CTick);
        ylabel(hC,'Approximation percentage error','fontsize',FS)
        line('parent',hC,'xdata',[-5,5],'ydata',log10(0.1  )*[1,1],'linewidth',LW,'color','k')
        line('parent',hC,'xdata',[-5,5],'ydata',log10(0.01 )*[1,1],'linewidth',LW,'color','k')
        line('parent',hC,'xdata',[-5,5],'ydata',log10(0.001)*[1,1],'linewidth',LW,'color','k')
        set(gca,'fontsize',FS)
    else
        xlabel('XLAB','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(0:3), 'xticklabel', {'X1','X2','X3','X4'})
        ylabel('YLAB','fontsize',FS)
        set(gca,'ylim',[0,1],'ytick',0:0.1:1, 'yticklabel', {'Y01','Y02','Y03','Y04','Y05','Y06','Y07','Y08','Y09','Y10','Y11'})
        caxis(log10(CTick([1,end])))
        colormap(cbrewer('seq','YlOrRd',50))
        hC = colorbar('ytick',log10(CTick),'yticklabel',{'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11','C12','C13'});
        ylabel(hC,'CLAB','fontsize',FS)
        line('parent',hC,'xdata',[-5,5],'ydata',log10(0.1  )*[1,1],'linewidth',LW,'color','k')
        line('parent',hC,'xdata',[-5,5],'ydata',log10(0.01 )*[1,1],'linewidth',LW,'color','k')
        line('parent',hC,'xdata',[-5,5],'ydata',log10(0.001)*[1,1],'linewidth',LW,'color','k')
        set(gca,'fontsize',FS)
    end
    axis square
    
    % Tidy
    clear CTick hC
end


% Threshold Lewis values - Liquid phase only
if 1 && SaltEffects
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Graph options
    CTick = [30,50,100,190,300,500,1000,2000,3000,5000,10000];
    
    % Create graph
    hold on
    contourf(MushyDiffVec,LqdusInfVec,log10(ThresholdLewis)',50,'edgecolor','none')
    contour( MushyDiffVec,LqdusInfVec,log10(ThresholdLewis)',log10([1000,170,100]),'edgecolor','k','linewidth',LW)
    patch(MushyDiffVec([1,1,end,end]),[0,0.05,0.05,0],0.5*[1,1,1])
    patch(MushyDiffVec([1,1,end,end]),[1,LqdusInfVec(end),LqdusInfVec(end),1],0.5*[1,1,1])
    hold off
    if 0
        xlabel('Mushy diffusivity factor,  m','fontsize',FS)
        set(gca,'xscale','log','xtick',10.^(0:3))
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
        set(gca,'xscale','log','xtick',10.^(0:3), 'xticklabel', {'X1','X2','X3','X4'})
        ylabel('YLAB','fontsize',FS)
        set(gca,'ylim',[0,1],'ytick',0:0.1:1, 'yticklabel', {'Y01','Y02','Y03','Y04','Y05','Y06','Y07','Y08','Y09','Y10','Y11'})
        caxis(log10(CTick([1,end])))
        colormap(cbrewer('seq','YlOrRd',50))
        hC = colorbar('ytick',log10(CTick),'yticklabel',CTick,'fontsize',FS,'yticklabel',{'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11'});
        ylabel(hC,'CLAB','fontsize',FS)
        line('parent',hC,'xdata',[-5,5],'ydata',log10( 100*[1,1]),'color','k','linewidth',LW)
        line('parent',hC,'xdata',[-5,5],'ydata',log10( 190*[1,1]),'color','k','linewidth',LW)
        line('parent',hC,'xdata',[-5,5],'ydata',log10(1000*[1,1]),'color','k','linewidth',LW)
        set(gca,'fontsize',FS)
        box on
    end
    
    % Tidy
    clear hC CTick
end

% Tidy
clear FC FS LW

