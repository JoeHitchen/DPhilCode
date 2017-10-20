function DiffMatt = DiffMatrix(NGrid,Depth)
% This function creates the Chebyshev differentiation matrix for NGrid+1
% points over a domain of size 'Depth'.
%
% The function implements formulas equivalent to transforming to modal
% space, performing a differentiation in modal space and transforming back
% to real space.
%
% Inputs:  NGrid - The number of grid points (Produces NGrid+1 points).
%          Depth - The domain depth.
%
% Outputs: D - The differentiation matrix.
%
% (17/04/15)
    
    % Prepare
    O = ones(NGrid+1,1);
    c = [2,ones(1,NGrid-1),2].*(-1).^(0:NGrid);
    X = cos((NGrid:-1:0)*pi/NGrid);
    
    % Non-diagonal elements
    DiffMatt = ((1./c)'*c)./(O*X-X'*O');
    
    % Diagonal elements
    DiffMatt(1:NGrid+2:(NGrid+1)^2) = O;
    DiffMatt(1:NGrid+2:(NGrid+1)^2) = DiffMatt(1:NGrid+2:(NGrid+1)^2) - sum(DiffMatt);
    
    % Scale by depth
    DiffMatt = DiffMatt*2/Depth;
    
    % Transpose (for some reason)
    DiffMatt = DiffMatt';
end
