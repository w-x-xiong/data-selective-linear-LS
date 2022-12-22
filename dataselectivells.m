function [x_tilde] = dataselectivells(Tx, Rx, Rg, n)
% DATA SELECTIVE LLS method
% by Wenxin Xiong
% Please refer to our paper
% ``Data-Selective Least Squares Methods for Elliptic Localization
% With NLOS Mitigation'' for more details.
% Rg: range sum of arrival measurement matrix (M*L)
% Tx: matrix holding the positions of M transmitters
% Rx: matrix holding the positions of L receivers
% n: user-specified number of TSOA measurements being used

[M, L] = size(Rg);

x_tilde_n = LLS(Tx, Rx, Rg);

x_tilde = x_tilde_n;

delta = 0;

for l = 1:L
    for m = 1:M
        delta = delta + (norm(x_tilde_n - Tx(:,m)) + norm(x_tilde_n - Rx(:,l)) - Rg(m,l))^2/(M*L);
    end
end

idx = zeros(2,M*L);
for m = 1:M
    for l = 1:L
        idx(1,M*(l-1)+m) = m;
        idx(2,M*(l-1)+m) = l;
    end
end

v = 1:M*L;
CMLn = nchoosek(v,n);


for i = 1:((factorial(M*L))/(factorial(n)*factorial(M*L-n)))
    
    Rg_new = zeros(M,L);
    
    CMLn_i = CMLn(i,:);
    
    for kk = 1:n
        Rg_new(idx(1, CMLn_i(kk)),idx(2, CMLn_i(kk))) = Rg(idx(1, CMLn_i(kk)),idx(2, CMLn_i(kk)));
    end
    
    A = zeros(M*L, L+2);

    q = zeros(M*L,1);
    
    for m = 1:M
        for l = 1:L
            A(M*(l-1)+m,:) = 2*[Rx(1,l)-Tx(1,m), Rx(2,l)-Tx(2,m), zeros(1,l-1), Rg(m,l), zeros(1,L-l)];
        end
    end

    for l = 1:L
        for m = 1:M
            q(M*(l-1)+m) = Rg(m,l)^2 + Rx(1,l)^2 + Rx(2,l)^2 - Tx(1,m)^2 - Tx(2,m)^2;
        end
    end
    
    A_new = [];
    
    q_new = [];
    
    for l = 1:L
        for m = 1:M
            if (Rg_new(m,l)~=0)
                A_new = [A_new;A(M*(l-1)+m,:)];
                q_new = [q_new;q(M*(l-1)+m)];
            end
        end   
    end
    
    Colnon0ele = any(A_new);
    
    A_new2 = [];
    
    for iii = 1:length(Colnon0ele)
        if (Colnon0ele(iii)~=0)
            A_new2 = [A_new2,A_new(:,iii)];
        end
    end
    
    y_tilde_tl = (pinv(A_new2'*A_new2))*A_new2'*q_new;

    x_tilde_tl = y_tilde_tl(1:2);
    
    f_msr = 0;
    
    for l = 1:L
        for m = 1:M
            if (Rg_new(m,l)~=0)
                f_msr = f_msr + (norm(x_tilde_tl - Tx(:,m)) + norm(x_tilde_tl - Rx(:,l)) - Rg(m,l))^2/n;
            end
        end
    end
    
    if f_msr < delta
        delta = f_msr;
        x_tilde = x_tilde_tl;
    end
    
end


end
