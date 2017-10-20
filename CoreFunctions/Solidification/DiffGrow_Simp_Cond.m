function InterfacePosition = DiffGrow_Simp_Cond(MushDiff,LqdusInf,varargin) 
% This function calculates the growth rate for the simplified model where
% the liquid fraction is large and the upper boundary is perfectly cooled.
% It is not capable of calculating the field profiles and cannot calculate
% ice depths at physics times.
%
% Inputs: MushDiff - Mushy diffusivity factor
%         LqdusInf - Initial liquidus temperature
% (Opt.)  SaltDiff - Salt diffusivity factor
%                      (Default value = Infinity)
% (Opt.)  IG       - Initial guess
%                      (Default value = 0.01)
%
% Output: InterfacePosition - Self similar interface position.
%
% (02/03/16)

    % Process input
    AddSalt = 0;
    if nargin > 2 && ~isinf(varargin{1})
        AddSalt = 1;
        SaltDiff = varargin{1};
    end
    
    % Specifiy initial guess
    if nargin > 3 && ~isempty(varargin{2})
        IG = varargin{2};
    else
        IG = 0.01;
    end
    
    % Set interface condition function
    IntCondFcn_NS = @(Pos) LqdusInf*erfcx(Pos/2) ...
        - (1-LqdusInf)*erf(MushDiff*Pos/2)*exp((MushDiff*Pos/2)^2)/MushDiff;
    
    if AddSalt
        IntCondFcn = @(Pos) IntCondFcn_NS(Pos) - erfcx(SaltDiff*Pos/2)/SaltDiff;
    else
        IntCondFcn = @(Pos) IntCondFcn_NS(Pos);
    end
    
    % Find interface
    [InterfacePosition, ~, EF] = fzero(IntCondFcn,IG,optimset('Display','off'));
    if EF ~= 1
        InterfacePosition = NaN;
    end
end

