% Copyright May. 25, 2019, Dr.WEN You-Wei
% email: wenyouwei@gmail.com

classdef ClassBlurMatrix
    properties
        psf
        eigBM
        boundarycond  = 'cir'
        tgv        = 0
    end
    
    methods
        function ob = ClassBlurMatrix(psf, imsize, flag)
            if nargin < 3, flag = 'cir';   end
            ob.boundarycond = flag;
            if length(imsize)>2 && imsize(3)>1 && size(psf,3) == 1
                psf  = repmat(psf,1,1,imsize(3));
            end
            ob.psf   = psf;
            ob.eigBM = zeros(imsize);
            
            s        = fix(([size(psf,1) size(psf,2)]+1)/2);
            switch ob.boundarycond
                case {'cir','circular'}
                    ob.eigBM(1:size(psf,1), 1:size(psf,2),:) = psf;
                    ob.eigBM = circshift(ob.eigBM, 1-s);
                    ob.eigBM = fft2(ob.eigBM);
                case {'refl','symmetric'}
                    psfa = psf(end:-1:1,end:-1:1,:);
                    if norm(psfa(:)-psf(:),'fro')>1e-7
                        error('The blur kernel should be symetric!');
                    end
                    h1 = psf(s(1):end, s(2):end,:);
                    h2 = h1(2:end,:,:); h2(s(1),:,:) = 0;
                    h3 = h1 + h2;
                    h4 = h3(:,2:end,:); h4(:,s(2),:) = 0;
                    h5 = h4 + h3;
                    
                    ob.eigBM(1:s(1),1:s(2),:) = h5;
                    e1 = zeros(imsize); e1(1,1,:) = 1;
                    for k = 1:imsize(3)
                        ob.eigBM(:,:,k) = dct2(ob.eigBM(:,:,k))./dct2(e1(:,:,k));
                    end
                    %ob.ForTF         = @(x)dct2(x);
                    %ob.BackTF        = @(x)idct2(x);
                otherwise
                    error('Wrong. Flag should be cir or refl or imfilter');
            end
        end
        
        function ob = ctranspose(A) %% written by wenyouwei@gmail.com
            ob       = A;
            ob.eigBM = conj(A.eigBM);
        end
        function y = mldivide(A,x)
            if ~isa(A,'ClassBlurMatrix')
                error('A must be the class of ClassBlurMatrix');
            end
            B = A; B.eigBM = 1./A.eigBM;
            y = B * x; 
        end
        
        function y = mrdivide(x,A)
            B = A; B.eigBM = 1./A.eigBM;
            y = B * x; 
        end
        
        function ob = abs(A) %% written by wenyouwei@gmail.com
            ob       = A;
            ob.eigBM = abs(A.eigBM);
        end
        
        function ob = inv(A) %% written by wenyouwei@gmail.com
            ob       = A;
            ob.eigBM = 1./(A.eigBM);
        end
        
        function ob = plus(a,b)
            if isa(a,'ClassBlurMatrix')
                ob = a;
                if isa(b,'ClassBlurMatrix')
                    ob.eigBM = a.eigBM + b.eigBM;
                else
                    ob.eigBM = a.eigBM + b;
                end
            else
                ob = b;
                ob.eigBM = a + b.eigBM;
            end
        end
        function ob = minus(a,b)
            if isa(a,'ClassBlurMatrix')
                ob = a;
                if isa(b,'ClassBlurMatrix')
                    ob.eigBM = a.eigBM - b.eigBM;
                else
                    ob.eigBM = a.eigBM - b;
                end
            else
                ob       = b;
                ob.eigBM = a - b.eigBM;
            end
        end
        
        
        function y = mtimes(a,x)%% written by wenyouwei@gmail.com
            if ~isa(a, 'ClassBlurMatrix')
                if numel(a) ~= 1
                    error('Wrong, 1th input: a scalar or a class');
                else
                    y = x;
                    y.eigBM = a * x.eigBM;
                end
                return;
            end
            if isa(x,'ClassBlurMatrix')
                y = a;    y.eigBM = y.eigBM .* x.eigBM;
                return;
            end
            
            if numel(x) == 1
                y = a;    y.eigBM = y.eigBM  * x;
                return;
            end
            switch a.boundarycond
                case {'cir','circular'}
                    y = ifft2(a.eigBM  .* fft2(x));
                    y = real(y); 
                case {'refl','symmetric'}
                    y = x; 
                    for j = 1:size(x,3)
                        y(:,:,j) = idct2(a.eigBM(:,:,j) .* dct2(x(:,:,j)));
                    end
            end
        end
        
    end
    
end

