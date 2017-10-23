%% Stefan Curve Investigation
% This script investigates the shape of the Stefan curves for a range of
% initial liquidus temperatures and liquidus gradients.
%
% Some numerical issues remain with the routine.
%
% It creates figures 3.12 (a&b) and 3.13 of my pre-corrections thesis.
%
% (14/01/16)

clear
clc
addpath(genpath('../CoreFunctions'))


%% Setup Investigation

% Investigation 1
I1_Run          = true;
I1_LqdusInfVec  = 0.5;
I1_StefanVec    = [0,10.^(-2:0.1:3)];
I1_LqdusGradVec = 10.^(-3:0.2:2);

% Investigation 4
I2_Run          = true;
%I2_LqdusInfVec  = [0.01:0.02:0.17,0.2:0.05:0.85,0.9:0.02:0.98,0.99];
I2_LqdusInfVec  = [0.01:0.02:0.17,0.2:0.05:0.85,0.9:0.02:0.96];
I2_StefanVec    = [0,10.^(-2:0.2:2)];
I2_LqdusGradVec = 10.^(-2.0:0.5:2);

% Computational Domain
Comp.DomainDepth = 15;
Comp.NGrid       = [39,49];


%% Run investigation 1

% Initialise result containers
I1_GrowthRate = NaN(length(I1_LqdusGradVec),length(I1_StefanVec));

% Loop over initial liquidus temperatures
if I1_Run
    disp('Running investigation 1...')
    
    % Loop over liquidus gradients
    hTic = tic;
    parfor iLg = 1:length(I1_LqdusGradVec)
        
        % Initialise results container
        PAR_GrowthRate = NaN(1,length(I1_StefanVec));
        
        % Loop over Stefan numbers
        for iSt = 1:length(I1_StefanVec)
            
            % Set IG
            if iSt == 1
                IG  = [DiffGrow_Simp_Cond(1,I1_LqdusInfVec),1/sqrt(pi),I1_LqdusGradVec(iLg)/(I1_LqdusGradVec(iLg)+I1_LqdusInfVec)];
            else
                IG(3) = I1_LqdusGradVec(iLg)/(I1_LqdusGradVec(iLg)+I1_LqdusInfVec);
            end
            
            % Get results
            [PAR_GrowthRate(iSt),~,~,IG] = DiffGrow_Mushy_Cond_NS(I1_StefanVec(iSt),I1_LqdusGradVec(iLg),I1_LqdusInfVec,Comp,IG);
        end
        
        % Output
        I1_GrowthRate(iLg,:) = PAR_GrowthRate;
        toc(hTic)
    end
    
    % Output
    disp('Investigation 1 complete!')
    disp(' ')
end

% Tidy
clear hTic iLg iSt PAR_Growth IG


%% Run investigation 4

% Initialise result containers
I2_GrowthRate = NaN(length(I2_LqdusInfVec),length(I2_StefanVec),length(I2_LqdusGradVec));

% Loop over initial liquidus temperatures
if I2_Run
    
    % Text output
    disp('Running investigation 2...')
    
    % Loop over liquidus gradients
    for iLg = 1:length(I2_LqdusGradVec)
        
        % Text output
        disp(['Liquidus gradient ',num2str(iLg),' of ',num2str(length(I2_LqdusGradVec))])
        
        % Loop over liquidus gradients
        hTic = tic;
        parfor iLi = 1:length(I2_LqdusInfVec)
%            disp([I2_LqdusGradVec(iLg), I2_LqdusInfVec(iLi)])
            
            % Initialise results container
            PAR_GrowthRate = NaN(1,length(I2_StefanVec));
    
            % Loop over Stefan numbers
            for iSt = 1:length(I2_StefanVec)
            
                % Set IG
                if iSt == 1
                    IG  = [DiffGrow_Simp_Cond(1,I2_LqdusInfVec(iLi)),1/sqrt(pi),I2_LqdusGradVec(iLg)/(I2_LqdusGradVec(iLg)+I2_LqdusInfVec(iLi))];
                else
                    IG(3) = I2_LqdusGradVec(iLg)/(I2_LqdusGradVec(iLg)+I2_LqdusInfVec(iLi));
                end
                
                % Get results
                [PAR_GrowthRate(iSt),~,~,IG] = DiffGrow_Mushy_Cond_NS(I2_StefanVec(iSt),I2_LqdusGradVec(iLg),I2_LqdusInfVec(iLi),Comp,IG);
            end
            
            % Output
            I2_GrowthRate(iLi,:,iLg) = PAR_GrowthRate;
            toc(hTic)
        end
        disp(' ')
    end
    
    % Output
    disp('Investigation 2 complete!')
    disp(' ')
end

% Tidy
clear hTic iLg iSt PAR_Growth IG


%% Stefan curve shifts

% Set options
StefanShift_UpperGrowthLim = 0.98;

% Get Stefan curve shape
StefanCurve_Stefan = [0,10.^(-7:0.1:7)];
StefanCurve_Growth = NaN(1,length(StefanCurve_Stefan));

for iSt = 1:length(StefanCurve_Stefan)
    
    % Calculate simplified growth rate
    StefanCurve_Growth(iSt) = DiffGrow_Simp_Cond(sqrt(1+StefanCurve_Stefan(iSt)),0.5);
end
StefanCurve_Growth = StefanCurve_Growth/StefanCurve_Growth(1);

% Calculate shifts
if I2_Run
    
    % Create results container
    I2_StefanShift = NaN(length(I2_LqdusInfVec),length(I2_LqdusGradVec));
    
    % Initial liquidus temperature
    for iLi = 1:length(I2_LqdusInfVec)
        for iLg = 1:length(I2_LqdusGradVec)
            
            % Select data
            St_Data = I2_StefanVec;
            GR_Data = I2_GrowthRate(iLi,:,iLg)/I2_GrowthRate(iLi,1,iLg);
            St_Data = St_Data(~isnan(GR_Data));
            GR_Data = GR_Data(~isnan(GR_Data));
            [~,UniqueIdx] = unique(GR_Data);
            St_Data = St_Data(UniqueIdx);
            GR_Data = GR_Data(UniqueIdx);
            
            % Calculate shift
            ShiftFcn = @(Growth) interp1(StefanCurve_Growth,StefanCurve_Stefan,Growth,'pchip')./interp1(GR_Data,St_Data,Growth,'pchip');
            if GR_Data(1) < StefanShift_UpperGrowthLim
                I2_StefanShift(iLi,iLg) = integral(ShiftFcn,GR_Data(1),StefanShift_UpperGrowthLim)/(StefanShift_UpperGrowthLim-GR_Data(1));
            end
        end
    end
end

% Perform fit
LIData = I2_LqdusInfVec'*ones(1,length(I2_LqdusGradVec));
LGData = ones(length(I2_LqdusInfVec),1)*I2_LqdusGradVec;
ShData = I2_StefanShift(:);

LIData = LIData(~isnan(ShData));
LGData = LGData(~isnan(ShData));
ShData = ShData(~isnan(ShData));


if 0
    FitFcn = @(Var) Var(1)*(LIData.^(5/3)).*(LGData.^(-2/3));
    Fit0 = [1];
elseif 0
    FitFcn = @(Var) Var(1)*(LIData.^(1+Var(2))).*(LGData.^(-Var(2)));
    Fit0 = [1,2/3];
else
    FitFcn = @(Var) Var(1)*(LIData.^(Var(3)+Var(2))).*(LGData.^(-Var(2)));
    Fit0 = [1,2/3,3/4];
end
ErrFcn = @(Var) sum(RelErr(FitFcn(Var),ShData,1).^2);
Fit1 = fminsearch(ErrFcn,Fit0);

ShiftFit1 = FitFcn(Fit0);
ShiftFit2 = FitFcn(Fit1);

% Tidy
clear iSt iLi iLg St_Data GR_Data UniqueIdx ShiftFcn StefanShift_UpperGrowthLim


%% Output
FC = 0;
LW = 3;
FS = 18;
PSFrag = true;

% Self-similar growth rate (mushy) surface plot
if 0 && I1_Run
    
    % Settings
    NContours = 50;
    CTicks = 0:0.1:1;
    RCont = [0.2:0.2:0.8,0.95];
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    Data = I1_GrowthRate'./(I1_GrowthRate(:,1)*ones(size(I1_StefanVec)))';
    
    % Create graph
    hold on
    contourf(I1_LqdusGradVec,I1_StefanVec(2:end),Data(2:end,:),NContours,'edgecolor','none')
    contour(I1_LqdusGradVec,I1_StefanVec(2:end),Data(2:end,:),[0.2:0.2:0.8,0.95],'edgecolor','r','linewidth',LW)
    hold off
    ylabel('Stefan number,  S_t','fontsize',FS)
    set(gca,'yscale','log','ytick',10.^(-2:3))
    xlabel('Liquidus gradient,  \Lambda','fontsize',FS)
    set(gca,'xscale','log','xlim',I1_LqdusGradVec([1,end]),'xtick',10.^(-3:2))
    set(gca,'fontsize',FS)
    caxis(CTicks([1,end]))
    hC = colorbar('ylim',CTicks([1,end]),'ytick',CTicks,'yticklabel',CTicks,'fontsize',FS);
    ylabel(hC,'Mushy growth rate,  \lambda')
    for Cont = RCont
        line('parent', hC, 'xdata', [-5,5], 'ydata', Cont*[1,1], 'color', 'r', 'linewidth', LW)
    end
    box on
    colormap(cbrewer('seq','YlGnBu',2*NContours))
    
    % Tidy
    clear NContours CTicks Data hC
end


% Collapsed curves
if 1 && I2_Run
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create data
    iC = 0;
    for iLi = 1:length(I2_LqdusInfVec)
        for iLg = 1:length(I2_LqdusGradVec)
     
            % Check for bad data
            if ~isnan(I2_StefanShift(iLi,iLg))
                
                % Increment counter
                iC = iC+1;
                
                % Select data
                Idx = ~isnan(I2_GrowthRate(iLi,:,iLg));
                Idx = Idx & ( I2_StefanVec*I2_StefanShift(iLi,iLg) > I2_StefanVec(2) );
                
                % Find curves
                ScaledCurves(:,iC) = interp1(I2_StefanVec(Idx)*I2_StefanShift(iLi,iLg),I2_GrowthRate(iLi,Idx,iLg)/I2_GrowthRate(iLi,1,iLg),StefanCurve_Stefan,'pchip',NaN); %#ok<SAGROW>
            end
        end
    end
    MaxCurve = max(ScaledCurves,[],2);
    MinCurve = min(ScaledCurves,[],2);
    
    % Create graph
    hold all
%    plot(StefanCurve_Stefan,ScaledCurves)
    hP(2) = plot(StefanCurve_Stefan,MaxCurve,'color',[1,0.5,0],'linewidth',LW);
    hP(3) = plot(StefanCurve_Stefan,MinCurve,'color',[1,0.0,0],'linewidth',LW);
    hP(1) = plot(StefanCurve_Stefan,StefanCurve_Growth/StefanCurve_Growth(1),'k','linewidth',LW);
    hold off
    if ~PSFrag
        xlabel('Scaled Stefan number','fontsize',FS)
        ylabel('Scaled growth rate','fontsize',FS)
        set(gca,'xscale','log','xlim',10.^[-2,5],'xtick',10.^(-2:5))
        set(gca,'ylim',[0,1],'ytick',0:0.1:1)
        set(gca,'fontsize',FS)
        legend(hP,'Model curve','Upper bound','Lower bound')
    else
        title('1TITLE','fontsize',FS)
        xlabel('1XLAB','fontsize',FS)
        ylabel('1YLAB','fontsize',FS)
        set(gca,'xscale','log','xlim',10.^[-2,5],'xtick',10.^(-2:5),'xticklabel',{'1X1','1X2','1X3','1X4','1X5','1X6','1X7','1X8'})
        set(gca,'ylim',[0,1],'ytick',0:0.2:1,'yticklabel',{'1Y1','1Y2','1Y3','1Y4','1Y5','1Y6'})
        set(gca,'fontsize',FS)
        legend(hP,'1LEG111111','1LEG2','1LEG3')
    end
    axis square
    box on
end


% Stefan shifts (calculated)
if 1 && I2_Run
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Create graph
    plot(I2_LqdusInfVec,I2_StefanShift(:,1:2:7),'-','linewidth',LW)
    if ~PSFrag
        xlabel('Initial liquidus temperature,  \theta_{L\infty}','fontsize',FS)
        set(gca,'xtick',0:0.1:1)
        ylabel('Stefan number prefactor,  p','fontsize',FS)
        set(gca,'yscale','log','ylim',10.^[-3,2],'ytick',10.^(-3:2))
        legend('\Lambda = 0.01','\Lambda = 0.1','\Lambda = 1','\Lambda = 10','location','southeast')
    else
        title('2TITLE','fontsize',FS)
        xlabel('2XLAB','fontsize',FS)
        set(gca,'xtick',0:0.2:1,'xticklabel',{'2X1','2X2','2X3','2X4','2X5','2X6'})
        ylabel('2YLAB','fontsize',FS)
        set(gca,'yscale','log','ylim',10.^[-3,2],'ytick',10.^(-3:2),'yticklabel',{'2Y1','2Y2','2Y3','2Y4','2Y5','2Y6'})
        legend('2LEG111','2LEG2','2LEG3','2LEG4','location','southeast')
    end
    set(gca,'fontsize',FS)
    axis square
    box on
end


% Stefan shifts (comparison)
if 1 && I2_Run
    
    % Create window
    FC = FC+1;
    figure(FC)
    clf
    
    % Percentage error calculation
    disp(   sqrt(sum(sum(((ShiftFit1-I2_StefanShift(~isnan(I2_StefanShift)))./I2_StefanShift(~isnan(I2_StefanShift))).^2))/length(I2_StefanShift(~isnan(I2_StefanShift))))   )
    
    % Create graph
    hold all
    plot(10.^[-3,2],10.^[-3,2],'k--','linewidth',LW)
    plot(I2_StefanShift(~isnan(I2_StefanShift)),ShiftFit1,'x','linewidth',LW)
    hold off
    if ~PSFrag
        xlabel('Measured prefactor,  p','fontsize',FS)
        set(gca,'xscale','log','xlim',10.^[-3,2])
        ylabel('Modelled prefactor,  p^{mod}','fontsize',FS)
        set(gca,'yscale','log','ylim',10.^[-3,2])
    else
        xlabel('XLAB','fontsize',FS)
        set(gca,'xscale','log','xlim',10.^[-3,2],'xtick',10.^(-3:2),'xticklabel',{'X1','X2','X3','X4','X5','X6'})
        ylabel('YLAB','fontsize',FS)
        set(gca,'yscale','log','ylim',10.^[-3,2],'ytick',10.^(-3:2),'yticklabel',{'Y1','Y2','Y3','Y4','Y5','Y6'})
    end
    set(gca,'fontsize',FS)
    axis square
    box on
end

% Tidy
clear FC FS LW


