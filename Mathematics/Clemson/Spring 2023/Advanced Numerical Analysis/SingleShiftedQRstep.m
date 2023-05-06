function [Hnew,HVs] = SingleShiftedQRstep(H,Mu2)

[m2,n2] = size(Mu2);
if m2 ~= 2 || n2 ~= 2 || ~isreal(Mu2)
    error('Wrong input Mu2: it must be a 2x2 real matrix.');
end
[m,n] = size(H);
if m ~= n
    error('Input H is not a square matrix!');
end
if nnz(tril(H,-2)) > 0
    error('Input H is not upper Hessenberg!');
end
if ~isreal(H)
    error('H is not real!');
end

HVs = zeros(3,n-1);

s = Mu2(1,1)+Mu2(2,2);
t = Mu2(1,1)*Mu2(2,2)-Mu2(1,2)*Mu2(2,1);

x = H(1,1)^2+H(1,2)*H(2,1)-s*H(1,1)+t;
y = H(2,1)*(H(1,1)+H(2,2)-s);
z = H(2,1)*H(3,2);
for k = 0 : n-3
    if y == 0 && z == 0
        continue;
    end
    signx = sign(x);
    if signx == 0,   signx = 1;  end
    vk = signx*norm([x; y; z])*eye(3,1) + [x; y; z];
    vk = vk/norm(vk);   vk = vk/norm(vk);
    HVs(:,k+1) = vk;
    q = max([1 k]);
    H(k+1:k+3,q:n) = H(k+1:k+3,q:n) - (2*vk)*(vk'*H(k+1:k+3,q:n));
    r = min([k+4 n]);
    H(1:r,k+1:k+3) = H(1:r,k+1:k+3) - (H(1:r,k+1:k+3)*(2*vk))*vk';
    x = H(k+2,k+1);
    y = H(k+3,k+1);
    if k < n-3
        z = H(k+4,k+1);
    end
end
if y ~= 0
    signx = sign(x);
    if signx == 0,  signx = 1;  end
    vk = signx*norm([x; y])*eye(2,1) + [x; y];
    vk = vk/norm(vk);   vk = vk/norm(vk);
    HVs(1:2,n-1) = vk;
    H(n-1:n,n-2:n) = H(n-1:n,n-2:n)-(2*vk)*(vk'*H(n-1:n,n-2:n));
    H(1:n,n-1:n) = H(1:n,n-1:n)-(H(1:n,n-1:n)*(2*vk))*vk';
end
Hnew = triu(H,-1);
