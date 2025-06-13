function [U,S,V] = svdecon(X,isPar)
% Input:
% X : m x n matrix
%
% Output:
% X = U*S*V'
%
% Description:
% Does equivalent to svd(X,'econ') but faster
%
% Vipin Vijayan (2014)

if nargin < 2
   isPar = false; 
end

%X = bsxfun(@minus,X,mean(X,2));
[m,n] = size(X);

if  m <= n
    if isPar
        spmd
            C = mtimes(X,X');
            [U,D] = eig(C);
            Ures = spmdCat(U,2,1);
            Dres = spmdCat(D,2,1);
        end
        U = Ures{1};
        D = Dres{1};
    else
        C = X*X';
        [U,D] = eig(C);
    end
    clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');
    U = U(:,ix);    
    
    if nargout > 2
        V = X'*U;
        s = sqrt(d);
        V = bsxfun(@(x,c)x./c, V, s');
        S = diag(s);
    end
else
    if isPar
        spmd
            C = mtimes(X',X);
            [V,D] = eig(C);
            Vres = spmdCat(V,2,1);
            Dres = spmdCat(D,2,1);
        end
        V = Vres{1};
        D = Dres{1};
    else
        C = X'*X; 
        [V,D] = eig(C);
    end
    clear C;
    
    [d,ix] = sort(abs(diag(D)),'descend');
    V = V(:,ix);    
    
    U = X*V; % convert evecs from X'*X to X*X'. the evals are the same.
    %s = sqrt(sum(U.^2,1))';
    s = sqrt(d);
    U = bsxfun(@(x,c)x./c, U, s');
    S = diag(s);
end
