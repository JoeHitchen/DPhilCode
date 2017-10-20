function Grid = ChebGrid(NGrid,Depth)
% This function creates a Chebyshev spaced grid.
%
% Inputs: NGrid - Number of grid points (Produces NGrid+1 points)
%         Depth - Domain depth
%
% (20/02/15)

    % Calculate grid
    Grid = ( cos( (NGrid:-1:0) * pi/NGrid ) +1 )' / (2/Depth);
end

