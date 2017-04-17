function p = probM(m,wms,smin,smax,varargin)
%% probM
%
%   p = probM(m,smin,smax)
%
%   Finds the probability of measurements by margnializing out the samples
%   that could have happened.
%
%%

%% Defaults
estimator_default.type = 'BLS';

method_opts_default.type = 'quad';
method_opts_default.dx = 10;


%% Parse inputs
Parser = inputParser;

addRequired(Parser,'m')
addRequired(Parser,'wms')
addRequired(Parser,'smin')
addRequired(Parser,'smax')
addParameter(Parser,'method_opts',method_opts_default)

parse(Parser,m,wms,smin,smax,varargin{:})

m = Parser.Results.m;
wms = Parser.Results.wms;
smin = Parser.Results.smin;
smax = Parser.Results.smax;
method_opts = Parser.Results.method_opts;


N = size(m,2);
if length(wms) ~= N && length(wms) == 1
    wms = wms*ones(1,N);
elseif length(wms) ~=N
    error('Number of entries in wms must equal number of measurements (number of columns in m)')
end

%% 

switch method_opts.type
    case 'quad'
        
        % Create x-vector
        dx = method_opts.dx;
        x = smin:dx:smax;
        
        % Create Simpson's nodes
        l = length(x);
        h = (smax - smin)/l;
        w = ones(1,l);
        w(2:2:l-1) = 4;
        w(3:2:l-1) = 2;
        w = w*h/3;
        
        % Reshape measurements for processing
        M = permute(m,[2 3 1]);
        M = repmat(M,[1,1,1,l]);
        x = reshape(x,[1 1 1 l]);
        X = repmat(x,[size(M,1) 1 size(M,3) 1]);
        
        % Generate estimate
        w = reshape(w,[1 1 1 l]);
        w = repmat(w,[1 1 size(m,1) 1]);
        likelihood = likelihoods(M,X,wms);
        p = sum(w.*likelihood,4);
        p = permute(p,[3 2 1]);
        
%         if length(wms) == 1
%             
%             % Create x-vector
%             dx = method_opts.dx;
%             x = smin:dx:smax;
%             
%             % Create Simpson's nodes
%             l = length(x);
%             h = (smax - smin)/l;
%             w = ones(1,l);
%             w(2:2:l-1) = 4;
%             w(3:2:l-1) = 2;
%             w = w*h/3;
%             
%             % Reshape measurements for processing
%             M = permute(m,[2 3 1]);
%             M = repmat(M,[1,1,1,l]);
%             x = reshape(x,[1 1 1 l]);
%             X = repmat(x,[size(M,1) 1 size(M,3) 1]);
%             
%             % Generate estimate
%             w = reshape(w,[1 1 1 l]);
%             w = repmat(w,[1 1 size(m,1) 1]);
%             likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%             p = sum(w.*likelihood,4);
%             p = permute(p,[3 2 1]);
%             
%         else
%             
%             if N == 1
%                 % Create x-vector
%                 dx = method_opts.dx;
%                 x = smin:dx:smax;
%                 
%                 % Create Simpson's nodes
%                 l = length(x);
%                 h = (smax - smin)/l;
%                 w = ones(1,l);
%                 w(2:2:l-1) = 4;
%                 w(3:2:l-1) = 2;
%                 w = w*h/3;
%                 
%                 % Reshape measurements for processing
%                 M = permute(m,[2 3 1]);
%                 M = repmat(M,[1,1,1,l]);
%                 x = reshape(x,[1 1 1 l]);
%                 X = repmat(x,[size(M,1) 1 size(M,3) 1]);
%                 
%                 % Generate estimate
%                 w = reshape(w,[1 1 1 l]);
%                 w = repmat(w,[1 1 size(m,1) 1]);
%                 likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-M).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
%                 p = sum(w.*likelihood,4);
%                 p = permute(p,[3 2 1]);
%                 
%             elseif N == 2
%                 
%                 % Create x-vector
%                 dx = method_opts.dx;
%                 x = smin:dx:smax;
%                 
%                 % Create Simpson'€™s nodes
%                 l = length(x);
%                 h = (smax - smin)/l;
%                 w = ones(1,l);
%                 w(2:2:l-1) = 4;
%                 w(3:2:l-1) = 2;
%                 w = w*h/3;
%                 
%                 % Reshape measurements for processing
%                 M = permute(m,[2 3 1]);
%                 M = repmat(M,[1,1,1,l]);
%                 x = reshape(x,[1 1 1 l]);
%                 X = repmat(x,[size(M,1) 1 size(M,3) 1]);
%                 
%                 % Generate estimate
%                 w = reshape(w,[1 1 1 l]);
%                 w = repmat(w,[1 1 size(m,1) 1]);
%                 l1 = ( (1./sqrt(2*pi)/wm1/X(1,:,:,:)) .* ...
%                     exp( -(X(1,:,:,:)-M(1,:,:,:)).^2 ./...
%                     (2*wm1.^2.*X(1,:,:,:).^2) ) );
%                 l2 = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .* ...
%                     exp( -(X(1,:,:,:)-M(2,:,:,:)).^2 ./...
%                     (2*wm.^2.*X(1,:,:,:).^2) ) );
%                 likelihood = l1.*l2;
%                 p = sum(w.*likelihood,4);
%                 p = permute(p,[3 2 1]);
%                 
%             else
%                 error('N > 2 not supported for noise model with differing Weber fractions!')
%             end
%         end
        
    otherwise
        error(['Integration method ' method_opts.type ' not supported!'])
end


%% Functions

%% Likelihood function
function l = likelihoods(M,X,wms)

ls = nan(size(X));
for wmi = 1:length(wms)
    wm = wms(wmi);
    ls(wmi,:,:,:) = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .* ...
        exp( -(X(1,:,:,:)-M(wmi,:,:,:)).^2 ./...
        (2*wm.^2.*X(1,:,:,:).^2) ) );
end
l = prod(ls,1);