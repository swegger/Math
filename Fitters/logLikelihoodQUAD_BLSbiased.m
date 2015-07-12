function logL = logLikelihoodQUAD_BLSbiased( WM,WY,B,N,x,y,xmin,xmax,dx)
%% ParameterLikelihoodFunction
%
%   logL = logLiklihoodQUAD_BLSbiased(wm,wy,b,N,x,y,xmin,xmax,dx)
%
%   Determines the log likelihood of the parameters wm, wy and b given data
%   in x and y
%
%%

for ii = 1:length(WM);
    for jj = 1:length(WY)
        for kk = 1:length(B);
            wm = WM(ii);
            wy = WY(jj);
            b = B(kk);
            disp([wm wy b])
            
            for i = 1:length(N)
                n = N{i};
                if xmin-5*wm*xmin < 0
                    m = 0:dx:xmax+5*wm*xmax;
                else
                    m = xmin-5*wm*xmin:dx:xmax+5*wm*xmax;
                end
                
                if n*length(m)^n > 4000000
                    error('Surpasing reasonable memory limits; suggest increasing dx or decreasing N')
                end
                
                % Set up Simpson's nodes
                l = length(m);
                w = ones(1,l);
                h = (m(end)-m(1))/l;
                w(2:2:l-1) = 4;
                w(3:2:l-1) = 2;
                w = w*h/3;
                
                W = w(:);
                for j = 2:n
                    W = W*w;
                    W = W(:);
                end
                
                Mtemp = cell(1,n);
                [Mtemp{:}] = ndgrid(m);
                M = zeros(l^n,n);
                for j = 1:n
                    M(:,j) = [Mtemp{j}(:)];
                end
                
                method_opts.type = 'quad';
                method_opts.dx = dx;
                fBLS = ScalarBayesEstimators(M,wm,xmin,xmax,'method',method_opts);
                X = repmat(x{i}',numel(fBLS),1);
                Y = repmat(y{i}',numel(fBLS),1);
                fBLS = repmat(fBLS,1,size(X,2));
                
                p_y_take_fBLS = (1./sqrt(2.*pi.*wy.^2.*fBLS.^2)) .* exp( -(Y - (fBLS+b)).^2./(2.*wy.^2.*fBLS.^2) );
                p_m_take_x = (1./sqrt(2.*pi.*wm.^2.*X.^2)).^n .* exp( -sum((repmat(permute(M,[1 3 2]),[1 size(X,2) 1])-repmat(X,[1 1 n])).^2,3)./(2.*wm.^2.*X.^2) );
                integrand = p_y_take_fBLS.*p_m_take_x;
                
                likelihood = W'*integrand;
                
                logLi(i) = -sum(log(likelihood));
                
                clear M
            end
            logL(ii,jj,kk) = sum(logLi);
            
        end
    end
end

end

