%% Initial Liquidus Investigation (Conducting, without salt)
% This script investigates the effects of the initial liquidus temperature
% against the liquidus gradient for St = 15.8.
%
% (14/01/16)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Setup Investigation

% Parameter space
LqdusInfVec  = 0.1:0.1:0.9;
StefanVec    = [0,10.^(-1:0.1:1.2)];
LqdusGradVec = 10.^(-3:0.2:2);

% Computational Domain
Comp.DomainDepth = 15;
Comp.NGrid       = [39,49];


%% Run investigation

% Initialise result containers
GrowthRate = NaN(length(LqdusGradVec),length(LqdusInfVec));
Grid       = NaN(Comp.NGrid(1)+Comp.NGrid(2)+1,length(LqdusGradVec),length(LqdusInfVec));
LqdFrac    = NaN(Comp.NGrid(1)+Comp.NGrid(2)+1,length(LqdusGradVec),length(LqdusInfVec));

Param.LqdusGrad = 0;
Param.Stefan    = 0;
Param.LqdusInf  = 0;

% Loop over initial liquidus temperatures
for iLi = 1:length(LqdusInfVec)
    disp(['Initial liquidus temperature ',num2str(iLi),' of ',num2str(length(LqdusInfVec))])
    
    % Loop over liquidus gradients
    hTic = tic;
    parfor iLg = 1:length(LqdusGradVec)
        
        % Loop over Stefan numbers
        for iSt = 1:length(StefanVec)
            
            % Set IG
            if iSt == 1
                IG  = [1,0.5,0];
            else
                IG(3) = LqdusGradVec(iLg)/(LqdusGradVec(iLg)+LqdusInfVec(iLi)); %#ok<PFBNS>
            end
    
            % Get results
            [GrowthRate(iLg,iLi),tmpProfM,tmpGridM,IG] = DiffGrow_Mushy_Cond_NS(StefanVec(iSt),LqdusGradVec(iLg),LqdusInfVec(iLi),Comp,IG);
            
            % Get profiles for heat capacity investigation
            if size(tmpProfM,1) == length(tmpGridM)
                Grid(   :,iLg,iLi) = tmpGridM;
                LqdFrac(:,iLg,iLi) = tmpProfM(:,5);
            end
        end
        toc(hTic)
    end
    disp(' ')
end

% Tidy
clear iLi iLg iSt hTic tmpGridS tmpProfS tmpGridM tmpProfM IG



%% Output
FC = 0;
LW = 3;
FS = 18;
PSFrag = true;


% Self-similar growth rate (mushy) surface plot
if 1
    
    % Settings
    NContours = 20;
    CTicks = [0.1,0.2,0.3,0.5,1,2];
    LGIndex = 1:length(LqdusGradVec)-5;
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    hold on
    contourf(LqdusGradVec(LGIndex),LqdusInfVec,log10(GrowthRate(LGIndex,:))',NContours,'edgecolor','none')
    contour(LqdusGradVec(LGIndex),LqdusInfVec,log10(GrowthRate(LGIndex,:))',log10(CTicks),'edgecolor','k','linewidth',LW)
    plot(0.3*(0.1:0.01:0.9),(0.1:0.01:0.9),'r','linewidth',LW)
    hold off
    if ~PSFrag
        ylabel('Initial Liquidus Temperature,  \theta_{L\infty}','fontsize',FS)
        set(gca,'ylim',[0.1,0.901],'ytick',LqdusInfVec)
        xlabel('Liquidus gradient,  \Lambda','fontsize',FS)
        set(gca,'xscale','log','xlim',LqdusGradVec(LGIndex([1,end])),'xtick',10.^(-3:2))
        set(gca,'fontsize',FS)
        caxis(log10(CTicks([1,end])))
        hC = colorbar('ylim',log10(CTicks([1,end])),'ytick',log10(CTicks),'yticklabel',CTicks,'fontsize',FS);
        ylabel(hC,'Mushy growth rate,  \lambda')
    else
        ylabel('YLAB','fontsize',FS)
        set(gca,'ylim',[0.1,0.901],'ytick',LqdusInfVec,'yticklabel',{'Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9'})
        xlabel('XLAB','fontsize',FS)
        set(gca,'xscale','log','xlim',LqdusGradVec(LGIndex([1,end])),'xtick',10.^(-3:2),'xticklabel',{'X1','X2','X3','X4','X5'})
        set(gca,'fontsize',FS)
        caxis(log10(CTicks([1,end])))
        hC = colorbar('ylim',log10(CTicks([1,end])),'ytick',log10(CTicks),'yticklabel',{'C1','C2','C3','C4','C5','C6'},'fontsize',FS);
        ylabel(hC,'CLAB')
        for iC = CTicks
            line('parent',hC,'xdata',[-5,5],'ydata',log10(iC)*[1,1],'color','k','linewidth',LW)
        end
    end
    box on
    axis square
    colormap(cbrewer('seq','YlGnBu',2*NContours))
    
    % Tidy
    clear NContours CTicks LGIndex hC
end


% Liquid fraction surface
if 1
    
    % Set Stefan number sample
    CTicks = 0:0.2:1;
    LISamples = [2,5,9];
    DepthTick{1} = [0,0.3001];
    DepthTick{2} = 0:0.05:0.6;
    DepthTick{3} = 0:0.1:1;
    
    for iLi = 1:length(LISamples)
        
        % Create window
        FC = FC+1;
        figure(FC)
        clf
        
        % Create graph
        hold on
        contourf(ones(Comp.NGrid(1)+1,1)*LqdusGradVec(6:end-10),Grid(1:Comp.NGrid(1)+1,6:end-10,LISamples(iLi)),LqdFrac(1:Comp.NGrid(1)+1,6:end-10,LISamples(iLi)),50,'edgecolor','none')
        contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec(6:end-10),Grid(1:Comp.NGrid(1)+1,6:end-10,LISamples(iLi)),LqdFrac(1:Comp.NGrid(1)+1,6:end-10,LISamples(iLi)),0.1:0.1:0.9,'edgecolor','k','linewidth',1)
        contour( ones(Comp.NGrid(1)+1,1)*LqdusGradVec(6:end-10),Grid(1:Comp.NGrid(1)+1,6:end-10,LISamples(iLi)),LqdFrac(1:Comp.NGrid(1)+1,6:end-10,LISamples(iLi)),0.5,'edgecolor','r','linewidth',LW)
        plot(LqdusGradVec,GrowthRate(:,LISamples(iLi)),'k-','linewidth',LW)
        hold off
        if ~PSFrag
            ylabel('Depth,  \eta','fontsize',FS)
            set(gca,'ydir','reverse')
            xlabel('Liquidus gradient,  \Lambda','fontsize',FS)
            set(gca,'xscale','log','xtick',10.^(-3:1))
            set(gca,'ylim',[min(DepthTick{iLi}),max(DepthTick{iLi})],'ytick',DepthTick{iLi})
            set(gca,'fontsize',FS)
            caxis(CTicks([1,end]))
            hC = colorbar('ylim',CTicks([1,end]),'ydir','normal','ytick',CTicks,'yticklabel',CTicks,'fontsize',FS);
            ylabel(hC,'Liquid fraction,  \chi')
            box on
            axis square
            colormap(cbrewer('seq','Blues',50))
            line('parent',hC,'xdata',[-5,5],'ydata',[0.5,0.5],'color','r','linewidth',LW)
        else
            title([num2str(iLi),'TITLE'],'fontsize',FS)
            ylabel('YLAB','fontsize',FS)
            set(gca,'ydir','reverse')
            xlabel('XLAB','fontsize',FS)
            set(gca,'xscale','log','xtick',10.^(-2:1),'xticklabel',{'X1','X2','X3'})
            set(gca,'ylim',[min(DepthTick{iLi}),max(DepthTick{iLi})],'ytick',0:0.1:1,'yticklabel',{'Y01','Y02','Y03','Y04','Y05','Y06','Y07','Y08','Y09','Y10','Y11'})
            set(gca,'fontsize',FS)
            colormap(cbrewer('seq','Blues',50))
            caxis(CTicks([1,end]))
            hC = colorbar('ylim',CTicks([1,end]),'ydir','reverse','ytick',CTicks,'yticklabel',{'C1','C2','C3','C4','C5','C6'},'fontsize',FS);
            ylabel(hC,'CLAB')
            box on
            axis square
            for iC = 0.1:0.1:0.9
                line('parent',hC,'xdata',[-5,5],'ydata',iC*[1,1],'color','k','linewidth',1)
            end
            line('parent',hC,'xdata',[-5,5],'ydata',[0.5,0.5],'color','r','linewidth',LW)
        end
    end
    
    % Tidy
    clear LISamples CTicks hC DepthTick iLi
end


% Tidy
clear FC FS LW