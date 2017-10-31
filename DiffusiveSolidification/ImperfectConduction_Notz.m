%% Notz Comparison
% 
% The script calculates the evolution of a single sea ice setup for
% conditions appropriate for comparison to the field work of Notz & Worster
% (2008).
%
% It creates figure 3.18 of my pre-corrections thesis.
% It also creates the plots in box 4 and 6 of my EGU poster.
%
% (07/07/16)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Parameters and settings

% Physical parameters
Set.DepthScale   = 10;        % Meters
Set.CoolingCoeff = 6.3;       % Watts per Kelvin per meter^2
Set.ThermCond    = 0.523;     % Watts per Kelvin per meter
Set.ThermDiff    = 1.2*10^-7; % Meters^2 per second
Set.LatentHeat   = 79;        % Kelvin
Set.LiquidusGrad = 0.085;     % Kelvin per ppt
Set.Salinity     = 35;        % ppt
Set.AirTemp      = -30;       % Kelvin
Set.WaterTemp    = -1;        % Kelvin

% Derived parameters
Set.TimeScale = Set.DepthScale^2/Set.ThermDiff/(60^2);                   % Hours
Set.TempScale = Set.WaterTemp-Set.AirTemp;                               % Kelvin
Set.Biot      = Set.CoolingCoeff*Set.DepthScale/Set.ThermCond;           % Dimensionless
Set.LqdusGrad = Set.LiquidusGrad*Set.Salinity/Set.TempScale;             % Dimensionless
Set.LqdusInf  = (-Set.LiquidusGrad*Set.Salinity-Set.AirTemp)/Set.TempScale; % Dimensionless
Set.Stefan    = Set.LatentHeat/Set.TempScale;                            % Dimensionless

% Simulation settings
PDE.StartBiot = 10^-0.5;
PDE.MushyDynamics = true;
PDE.Depth = 1;
if 1
    PDE.NGrid = 400;
    PDE.RelTol = 10^-4;
    PDE.AbsTol = 10^-6;
else
    PDE.NGrid = 300;
    PDE.RelTol = 10^-4;
    PDE.AbsTol = 10^-6;
end
PDE.tGrid = [0:0.02:144]/Set.TimeScale;
PDE.zGrid = PDE.Depth*(0:PDE.NGrid).^2/PDE.NGrid^2;

% Ratios
Set.Ratios.ThermCond = 4.24;
Set.Ratios.HeatCap   = 0.501;

%% Freezing 

% Calculate freezing Biot
InterfaceTempFcn = @(Biot) Set.LqdusInf-BiotCoolingProfile(0,1,Biot,1);
PDE.FreezingBiot = fzero(InterfaceTempFcn,1);

% Perform integration
tic
[PDE.IntPos,PDE.Temps] = DiffGrow_PDE_Routine(Set.Stefan,Set.LqdusGrad,Set.LqdusInf,Set.Biot,PDE,Set.Ratios);
toc

% Get liquid fraction
PDE.LqdFracs = Set.LqdusGrad./(Set.LqdusGrad+Set.LqdusInf-PDE.Temps);
PDE.LqdFracs(PDE.LqdFracs>1) = 1;
PDE.LqdFracs(PDE.LqdFracs<0) = 1;

% Tidy
clear InterfaceTempFcn


%% Output
FC = 0;
LW = 3;
FS = 18;
FW = 'normal';
MS = 12;
PSFrag = false;

% Ice profiles
if 1
    
    % Settings
    TimeAxis  = 0:24:144;
    TempAxis  = -20:10:-0;
    SFAxis    = 0.1:0.2:0.5;
    DepthAxis = 0:0.05:0.25;
    NotzTime  = [6,12,24,48,72,96,138];
    NotzDepth = [3.8,5.8,9.8,13,17,NaN,NaN]*0.01;
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    %%% Create subplot
    subplot(3,1,1)
    
    % Create graph
    [hAx,hL1,hL2] = plotyy(Set.TimeScale*PDE.tGrid,PDE.Temps(1,:)*Set.TempScale+Set.AirTemp,Set.TimeScale*PDE.tGrid,PDE.LqdFracs(1,:));
    set(hL1,'linewidth',LW)
    set(hL2,'linewidth',LW)
    if ~PSFrag
        xlabel(hAx(1),'Time,  hours','fontsize',FS,'fontweight',FW)
        set(hAx(1),'xlim',TimeAxis([1,end]),'xtick',TimeAxis)
        ylabel(hAx(1),'Surf. Temp.,  ^oC','fontsize',FS,'fontweight',FW)
        set(hAx(1),'ylim',TempAxis([1,end]),'ytick',TempAxis)
        set(hAx(1),'ycolor','k','fontsize',FS,'fontweight',FW)
        set(hAx(2),'xaxislocation','top','xlim',TimeAxis([1,end]),'xtick',TimeAxis,'xticklabel',[])
        ylabel(hAx(2),'Surf. Porosity','fontsize',FS,'fontweight',FW)
        set(hAx(2),'yaxislocation','right','ylim',SFAxis([1,end]),'ytick',SFAxis)
        set(hAx(2),'ycolor','k','fontsize',FS,'fontweight',FW)
        legend([hL1,hL2],'Surface Temperature','Surface Porosity')
    else
        xlabel(hAx(1),'1XLAB','fontsize',FS,'fontweight',FW)
        set(hAx(1),'xlim',TimeAxis([1,end]),'xtick',TimeAxis,'xticklabel',{'1X1','1X2','1X3','1X4','1X5','1X6','1X7'})
        ylabel(hAx(1),'1YLAB','fontsize',FS,'fontweight',FW)
        set(hAx(1),'ylim',TempAxis([1,end]),'ytick',TempAxis,'yticklabel',{'1Y1','1Y2','1Y3'})
        set(hAx(1),'ycolor','k','fontsize',FS,'fontweight',FW)
        set(hAx(2),'xaxislocation','top','xlim',TimeAxis([1,end]),'xtick',TimeAxis,'xticklabel',[])
        ylabel(hAx(2),'2YLAB','fontsize',FS,'fontweight',FW)
        set(hAx(2),'yaxislocation','right','ylim',SFAxis([1,end]),'ytick',SFAxis,'yticklabel',{'2Y1','2Y2','2Y3'})
        set(hAx(2),'ycolor','k','fontsize',FS,'fontweight',FW)
        legend([hL1,hL2],'LEG1111111111111111','LEG2')
    end
    box off
    
    %%% Create subplot 2
    subplot(3,1,2:3)
    
    % Create graph
    Ax3 = gca;
    hold all
    contourf(Set.TimeScale*PDE.tGrid,Set.DepthScale*PDE.zGrid,PDE.LqdFracs,0:0.1:1)
    contour( Set.TimeScale*PDE.tGrid,Set.DepthScale*PDE.zGrid,PDE.LqdFracs,0.5,'color',[1.0,0.5,0],'linewidth',LW)
    plot(    Set.TimeScale*PDE.tGrid,[0,PDE.IntPos(2:end)]*Set.DepthScale,'r-','linewidth',LW)
    plot(NotzTime,NotzDepth,'x','color',[0,0.5,0],'linewidth',LW,'markersize',MS)
    hold off
    if ~PSFrag
        set(Ax3,'xlim',TimeAxis([1,end]),'xtick',TimeAxis,'xaxislocation','top')
        ylabel('Depth','fontsize',FS,'fontweight',FW)
        set(Ax3,'ydir','reverse','ylim',DepthAxis([1,end]),'ytick',DepthAxis)
    
        caxis([0,1])
        colormap(cbrewer('seq','Blues',10))
        hC = colorbar('location','southoutside','xdir','reverse','xlim',[0,1],'xtick',0:0.1:1,'fontweight',FW);
        xlabel(hC,'Liquid fraction, \chi','fontsize',FS)
    else
        set(Ax3,'xlim',TimeAxis([1,end]),'xtick',TimeAxis,'xaxislocation','top','xticklabel',{'3X1','3X2','3X3','3X4','3X5','3X6','3X7'})
        ylabel('3YLAB','fontsize',FS,'fontweight',FW)
        set(Ax3,'ydir','reverse','ylim',DepthAxis([1,end]),'ytick',DepthAxis,'yticklabel',{'3Y1','3Y2','3Y3','3Y4','3Y5','3Y6'})
    
        caxis([0,1])
        colormap(cbrewer('seq','Blues',10))
        hC = colorbar('location','southoutside','xdir','reverse','xlim',[0,1],'xtick',0:0.1:1,'xticklabel',{'C01','C02','C03','C04','C05','C06','C07','C08','C09','C10','C11'},'fontweight',FW);
        xlabel(hC,'CLAB','fontsize',FS)
    end
    for iC = 0:0.1:1
        line('parent',hC,'xdata',iC*[1,1],'ydata',[-5,5],'color','k','linewidth',1)
    end
    line('parent',hC,'xdata',min(PDE.LqdFracs(1,PDE.tGrid*Set.TimeScale<TimeAxis(end)))*[1,1],'ydata',[-5,5],'color','k','linewidth',LW)
    line('parent',hC,'xdata',1/(1+Set.LqdusInf/Set.LqdusGrad)*[1,1],'ydata',[-5,5],'color',0.5*[1,1,1],'linewidth',LW)
    line('parent',hC,'xdata',0.5*[1,1],'ydata',[-5,5],'color',[1,0.5,0],'linewidth',LW)
    line('parent',hC,'xdata',0.995*[1,1],'ydata',[-5,5],'color','r','linewidth',LW)
    set(Ax3,'fontsize',FS,'fontweight',FW)
    
    Ax4 = axes('position',get(Ax3,'position'),'color','none');
    set(Ax4,'xlim',TimeAxis([1,end]),'xtick',TimeAxis,'xticklabel',[])
    set(Ax4,'ydir','reverse','ylim',DepthAxis([1,end]),'ytick',DepthAxis,'yticklabel',[])
    box on
    
    % Tidy
    clear hAx hL1 hL2 Ax3 Ax4 hC TimeAxis TempAxis SFAxis DepthAxis NotzTime NotzDepth
end


% Tidy
clear FC LW FS FW MS