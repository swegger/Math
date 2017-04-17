function p = probM(M,wm,smin,smax,varargin)
%% probM
%
%   p = probM(M,smin,smax)
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

addRequired(Parser,'M')
addRequired(Parser,'wm')
addRequired(Parser,'smin')
addRequired(Parser,'smax')
addParameter(Parser,'estimator',estimator_default)
addParameter(Parser,'method_opts',method_opts_default)
addParameter(Parser,'wm1',NaN)

parse(Parser,M,wm,smin,smax,varargin{:})

M = Parser.Results.M;
wm = Parser.Results.wm;
smin = Parser.Results.smin;
smax = Parser.Results.smax;
estimator = Parser.Results.estimator;
method_opts = Parser.Results.method_opts;
wm1 = Parser.Results.wm1;


N = size(M,2);

%% 

switch method_opts.type
    case 'quad'
        if isnan(wm1)
            
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
            m = permute(M,[2 3 1]);
            m = repmat(m,[1,1,1,l]);
            x = reshape(x,[1 1 1 l]);
            X = repmat(x,[size(m,1) 1 size(m,3) 1]);
            
            % Generate estimate
            w = reshape(w,[1 1 1 l]);
            w = repmat(w,[1 1 size(M,1) 1]);
            likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-m).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
            p = sum(w.*likelihood,4);
            p = permute(p,[3 2 1]);
            
        else
            
            if N == 1
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
                m = permute(M,[2 3 1]);
                m = repmat(m,[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(m,1) 1 size(m,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(M,1) 1]);
                likelihood = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)).^N .* exp( -(sum((X-m).^2,1))./(2*wm.^2.*X(1,:,:,:).^2) ) );
                p = sum(w.*likelihood,4);
                p = permute(p,[3 2 1]);
                
            elseif N == 2
                
                % Create x-vector
                dx = method_opts.dx;
                x = smin:dx:smax;
                
                % Create Simpson'€™s nodes
                l = length(x);
                h = (smax - smin)/l;
                w = ones(1,l);
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                % Reshape measurements for processing
                m = permute(M,[2 3 1]);
                m = repmat(m,[1,1,1,l]);
                x = reshape(x,[1 1 1 l]);
                X = repmat(x,[size(M,1) 1 size(m,3) 1]);
                
                % Generate estimate
                w = reshape(w,[1 1 1 l]);
                w = repmat(w,[1 1 size(m,1) 1]);
                l1 = ( (1./sqrt(2*pi)/wm1/X(1,:,:,:)) .* ...
                    exp( -(X(1,:,:,:)-m(1,:,:,:)).^2 ./...
                    (2*wm1.^2.*X(1,:,:,:).^2) ) );
                l2 = ( (1./sqrt(2*pi)/wm/X(1,:,:,:)) .* ...
                    exp( -(X(1,:,:,:)-m(2,:,:,:)).^2 ./...
                    (2*wm.^2.*X(1,:,:,:).^2) ) );
                likelihood = l1.*l2;
                p = sum(w.*likelihood,4);
                p = permute(p,[3 2 1]);
                
            else
                error('N > 2 not supported for noise model with differing Weber fractions!')
            end
        end
        
    otherwise
        error(['Integration method ' method_opts.type ' not supported!'])
end