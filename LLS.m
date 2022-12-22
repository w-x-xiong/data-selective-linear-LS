function [x_tilde] = LLS(Tx, Rx, Rg)
%LLS in the position handbook
%Rg=zeros(M,L);
%for m=1:M%%Txs
%    for l=1:L  %%Rxs
%        Rg(m,l) = rserr((l-1)*M+m);
%    end
%end
[M, L] = size(Rg);

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

y_tilde = (pinv(A'*A))*A'*q;

x_tilde = y_tilde(1:2);


end

