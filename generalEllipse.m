function [x,y] = generalEllipse(primaryR,secondaryR,varargin)
%% generalEllipse
%
%   [x,y] = generalEllipse(primaryR,secondaryR)
%
%   Generates an ellipse with radius primaryR in the x-direction and radius
%   secondaryR in the y direction with a default angular resolution of pi/6
%   centered at 0,0.
%
%   [x,y] = generalEllipse(...,'Center',[x0 y0])
%
%   Generates an ellipse as above but centered at x0,y0.
%
%   [x,y] = generalEllipse(...,'theta',theta)
%
%   Generates an ellipse as above but at the angles specified by the vector
%   theta.
%
%   [x,y] = generalEllipse(...,'phi',phi)
%
%   Generates an ellipse as above but rotated so that the primary axis is
%   phi radians from x.
%
%%

% Defaults
theta_default = 0:pi/6:2*pi;
phi_default = 0;
x0_default = 0;
y0_default = 0;
plt_default = 0;

% Parse inputs
Parser = inputParser;
addRequired(Parser,'primaryR');
addRequired(Parser,'secondaryR');
addParameter(Parser,'Center',[x0_default y0_default]);
addParameter(Parser,'theta',theta_default);
addParameter(Parser,'phi',phi_default);
addParameter(Parser,'Plot',plt_default);

parse(Parser,primaryR,secondaryR,varargin{:})

primaryR = Parser.Results.primaryR;
secondaryR = Parser.Results.secondaryR;
Center = Parser.Results.Center;
theta = Parser.Results.theta;
phi = Parser.Results.phi;
Plot = Parser.Results.Plot;

% Generate the ellipse, before rotation
x = primaryR*cos(theta);
y = secondaryR*sin(theta);

% Rotate to the desired angle
if phi ~= 0
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];      % Rotation matrix
    new = [x;y]'*R;
    x = new(:,1);
    y = new(:,2);
end

% Push to desired center
x = x+Center(1);
y = y+Center(2);

% Plot, if requested
if Plot
    plot(x,y,'k');
end