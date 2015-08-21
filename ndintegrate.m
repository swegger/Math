function [integraldx, errest] = ndintegrate(functionHandle,MinMax,varargin)
%% ndintegrate
%
%   Computes the approximate integral over a given volume.
%
%   integraldx = ndintegrate(functionHandle,MinMax)
%       Integrates the function "functionHandle" over the volume specified
%       by MinMax using Simpson's quadrature. MinMax should be arranged as:
%               MinMax(1,:) = [Min_dimension_1 Max_dimension_1]
%               MinMax(2,:) = [Min_dimension_2 Max_dimension_2]
%               ...
%               Minmax(n,:) = [Min_dimension_n Max_dimension_n];
%       The function should be arranged so that inputs are in columns
%       and the output is a vector with a length equal to the number of
%       rows in the input.
%
%   integraldx = ndintegrate(functionHandle,MinMax,'method',method)
%       Integrates the function "functionHandle" over the volume specified
%       according to the method specified by "method". Available methods:
%           'quad' (default) - uses the composite Simpson's quadrature
%           'trapz' - uses linear interpolation
%           'MonteCarlo' - uses Monte Carlo integration
%           'Gauss-Lobatto' - NOT YET SUPPORTED
%           'Gauss' - TODO - NOT YET SUPPORTED
%           'Newton-Cotes' - NOT YET SUPPORTED
%
%   integraldx =
%   ndintegrate(functionHandle,MinMax,'method','quad','options',options)
%       Integrates the function "functionHandle" over the volume specified
%       by MinMax using Simpson's quadrature. Step size of the integration
%       is specified by the vector options.dx, where options.dx(i) is the
%       percent of the distance between MinMax(i,1) and MinMax(i,2)
%
%   integraldx =
%   ndintegrate(functionHandle,MinMax,'method','MonteCarlo','options',options)
%       Calculates the integral according to Monte Carlo integration using
%       N samples as specified by options.N;
%
%   [integraldx, error] = ndintegrate(...)
%       Calculates an estimate of the error of the integration. For
%       'trapz', 'quad', ..., the estimate is an upper bound on the error.
%       For 'MonteCarlo', the value of the error represents the 95%
%       confidence bound (default), or confidence bound supplied by the
%       user (options.bound; TODO);
%
%   integraldx = ndintegrate(...,'ExtraVaribles',ExtraVariables)
%       Allows the user to supply extra variables for the user supplied
%       funciton specified by functionHandle. Note, the user must carefully
%       plan the function so that the output for each different value of
%       the extra variables will end up as a matrix, M, with size(M,1)
%       equal to size(X,1) and the columns of M representing the output of
%       the function with different values of the extra variables.
%
%
%   by Seth Egger
%       Update       Author       Note
%       2014/05/07   swe          Initial commit
%       2014/06/19   swe          fixed Monte Carlo random point generator
%
%%

% defaults
options.dx = 0.01;      % Sets step size for Simpson's quadrature or trapz methods to 1% of MinMax(2,:)-MinMax(1,:)

% Parse inputs
p = inputParser;

addRequired(p,'functionHandle');
addRequired(p,'MinMax');
addParameter(p,'method','quad');
addParameter(p,'options',options);
addParameter(p,'ExtraVariables',[]);
addParameter(p,'OpenMind',false)

parse(p,functionHandle,MinMax,varargin{:});

functionHandle = p.Results.functionHandle;
MinMax = p.Results.MinMax;
method = p.Results.method;
options = p.Results.options;
y = p.Results.ExtraVariables;
OpenMind = p.Results.OpenMind;


% Perform the integration
switch method
    case 'quad'
        % Number of dimensions
        N = size(MinMax,1);
        
        % Set up vectors to integrate over
        dx = options.dx(:).*(diff(MinMax,1,2));
        x = cell(1,N);
        lx = nan(1,N);
        for i = 1:N
            x{i} = MinMax(i,1):dx(i):MinMax(i,2);
            lx(i) = length(x{i});
        end
        
        if sum(lx.^N) > 4000000 && ~OpenMind
            error('Likely to stress memory limits; suggest increasing options.dx or trying the Monte Carlo integration method')
        end
        
        % Set up Simpson's nodes
        w = cell(1,N);
        for i = 1:N
            w{i} = ones(1,lx(i));
            h = (x{i}(end)-x{i}(1))/lx(i);
            w{i}(2:2:lx(i)-1) = 4;
            w{i}(3:2:lx(i)-1) = 2;
            w{i} = w{i}*h/3;
        end
        
        W = w{1}(:);
        for i = 2:N
            W = W*w{i};
            W = W(:);
        end
        
        Xtemp = cell(1,N);
        [Xtemp{:}] = ndgrid(x{:});
        X = nan(length(Xtemp{1}(:)),N);
        for i = 1:N
            X(:,i) = Xtemp{i}(:);
        end
        
        if isempty(y)
            integrand = functionHandle(X);
        else
            integrand = functionHandle(X,y);
        end
        
        % Perform the integral
        integraldx = W'*integrand;
        
        % TODO
        % error E <= max( h^5/90 * f''''(x) )
        errest = NaN;
        
    case 'trapz'
        % Number of dimensions
        N = size(MinMax,1);
        
        % Set up vectors to integrate over
        dx = options.dx(:).*(diff(MinMax,1,2));
        x = cell(1,N);
        lx = nan(1,N);
        for i = 1:N
            x{i} = MinMax(i,1):dx(i):MinMax(i,2);
            lx(i) = length(x{i});
        end
        
        if sum(lx.^N) > 4000000 && ~OpenMind
            error('Likely to stress memory limits; suggest increasing options.dx or trying the Monte Carlo integration method')
        end        
        
        % Set up weights for linear interpolation  % CHECK IF THIS IS RIGHT
        w = cell(1,N);
        for i = 1:N
            w{i} = ones(1,lx(i));
            h = (x{i}(end)-x{i}(1))/lx(i);
            w{i}(2:lx(i)-1) = 2;
            w{i} = w{i}*h/2;
        end
        
        W = w{1}(:);
        for i = 2:N
            W = W*w{i};
            W = W(:);
        end
        
        Xtemp = cell(1,N);
        [Xtemp{:}] = ndgrid(x{:});
        X = nan(length(Xtemp{1}(:)),N);
        for i = 1:N
            X(:,i) = Xtemp{i}(:);
        end
        
        if isempty(y)
            integrand = functionHandle(X);
        else
            integrand = functionHandle(X,y);
        end
        
        % Perform the integral
        integraldx = W'*integrand;
        
        % TODO
        % No estimation of error
        errest = NaN;
            
        
    case 'MonteCarlo_batch'
        % Read off options or use defaults if none exist
        if isfield(options,'N')
            N = options.N;      % Number of iterations to find
        else
            N = 1000;
        end
        if isfield(options,'bound')
            bound = options.bound;
        else
            bound = 0.95;
        end
        if isfield(options,'batch_sz')
            batch_sz = options.batch_sz;
        else
            batch_sz = 10000;
        end
        
        nd = size(MinMax,1);    % Number of dimensions overwhich to integrate
        
        % Calculate the volume
        V = prod(diff(MinMax,1,2));
        
        ntot = 0;
        sfx = 0;
        for i = 1:ceil(N/batch_sz)
            if ntot + batch_sz > N
                n = N - ntot;
            else
                n = batch_sz;
            end
            ntot = ntot + n;
            
            % Generate samples from the volume specified by MinMax;
            x = repmat(diff(MinMax,1,2)',n,1).*rand(n,nd) + repmat(MinMax(:,1)',n,1);
            
            % Evaluate the function at each set of inputs
            if isempty(y)
                output = functionHandle(x);
            else
                output = functionHandle(x,y);
            end
            
            % Estimate of the average value of f(x) over the volume
            sfx = sfx + sum(output);
        end
        
        % Estimate the average value of f(x) over the volume
        mfx = sfx/N;
        
        % Sample variance
        vfx = var(output);
        
        % Estimate the integral
        integraldx = V*mfx;
        
        % Calculate confidence interval
        estvar = V^2*(vfx)/N;       % Variance of the estimate
        if bound == 0.95
            m = 1.96;
        else
            m = findm(bound);
        end
        errest = m*sqrt(estvar);
        
    case 'MonteCarlo'
        % Read off options or use defaults if none exist
        if isfield(options,'N')
            N = options.N;      % Number of iterations to find
        else
            N = 1000;
        end
        if isfield(options,'bound')
            bound = options.bound;
        else
            bound = 0.95;
        end
        
        nd = size(MinMax,1);    % Number of dimensions overwhich to integrate
        
        % Calculate the volume
        V = prod(diff(MinMax,1,2));
        
        % Generate samples from the volume specified by MinMax;
        x = repmat(diff(MinMax,1,2)',N,1).*rand(N,nd) + repmat(MinMax(:,1)',N,1);
        
        % Evaluate the function at each set of inputs
        if isempty(y)
            output = functionHandle(x);
        else
            output = functionHandle(x,y);
        end
        
        % Estimate of the average value of f(x) over the volume
        mfx = mean(output);
        
        % Sample variance
        vfx = var(output);
        
        % Estimate the integral
        integraldx = V*mfx;
        
        % Calculate confidence interval
        estvar = V^2*(vfx)/N;       % Variance of the estimate
        if bound == 0.95
            m = 1.96;
        else
            m = findm(bound);
        end
        errest = m*sqrt(estvar);
        
        
        
    case 'Gauss-Lobatto'
        
        error('Not yet supported')   
        
    case 'Gauss'
        
        error('Not yet supported')    
        
    case 'Newton-Cotes'
        
        error('Not yet supported')    
end

%% Functions

% findm TODO: improve findm to exact solution if possible
function m = findm(bound)
f = @(m)(abs(bound - (normcdf(m)-normcdf(-m))));
m = fminsearch(f,1);
