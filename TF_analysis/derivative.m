%==========================================================================
%
% derivative  Numerical differentiation of data (i.e. arrays) over the 
% domain of the data or at specified points.
%
%   dy = derivative(x,y)
%   dy = derivative(x,y,x_star)
%
% See also diff, gradient, iderivative, derivest.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2021-08-27
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/documentation/Basic_Numerical_Calculus.pdf
%
% REFERENCES:
%   [1] http://www.ohiouniversityfaculty.com/youngt/IntNumMeth/lecture27.pdf
%   [2] https://en.wikipedia.org/wiki/Finite_difference_method
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - ((N+1)×1 or 1×(N+1) double) independent variable data
%   y       - ((N+1)×1 or 1×(N+1) double) dependent variable data
%   x_star  - (OPTIONAL) (p×1 or 1×p double) points at which to 
%             differentiate
%
% -------
% OUTPUT:
% -------
%   dy      - ((N+1)×1, 1×(N+1), p×1, or 1×p double) derivative of y = f(x)
%             w.r.t. x evaluated at:
%               --> all the points in x (cumulative differentiation)
%               --> all the points in x_star (point differentiation)
%
% -----
% NOTE:
% -----
%   --> If "x_star" is not input, then "dy" stores the derivative of 
%       y = f(x) with with respect to x at the points in "x".
%   --> N = number of data points (i.e. length of "y" and "x")
%   --> p = number of points to differentiate at (i.e. length of "x_star")
%
%==========================================================================
function [dy,x] = derivative(x,y,x_star)

    % number of subintervals
    N = length(x)-1;

    % preallocates vector to store derivative
    dy = zeros(size(x));

    % approximates derivative at lower bound using forward difference
    dy(1) = (y(2)-y(1))/(x(2)-x(1));

    % approximates derivative at upper bound using backward difference
    dy(N+1) = (y(N+1)-y(N))/(x(N+1)-x(N));

    % approximates derivatives at all other nodes using central differences
    for i = 2:N
        dy(i) = (y(i+1)-y(i-1))/(x(i+1)-x(i-1));
    end
    
    % approximates derivative at specified points via linear interpolation
    if nargin == 3
        dy = interp1(x,dy,x_star,'linear','extrap');
    end
        
end