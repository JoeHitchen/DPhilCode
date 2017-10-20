function [Error] = RelErr(Test,True,varargin)
% This function calculates the relative error of the data in 'Test' when
% compared to the true values in 'True'
%
% Inputs: Test    - The values to be tested for their accuracy.
%         True    - The true values used to be in the comparison.
% (Opt.)  AbsFlag - Return absolute value if 1. 
%
% Output: Error - The relative errors
%
% (01/04/15)

    % Calculate relative errors
    Error = (Test-True)./True;
    
    % Calculate absolute value
    if (nargin > 2) && varargin{1}
        Error = abs(Error);
    end
end

