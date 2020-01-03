function Proj=calcord(N)

Ix=.5*[0  1; 1 0];
A=ones(2^N);
%A=kron(Ix,eye(2));

n=2*N+1;

Iz=.5*[1  0; 0 -1];

Fz=Iz;
for i=2:N
  Fz=kron(Fz,eye(2))+kron(eye(length(Fz)),Iz);
end
R=expm(-j*Fz*2*pi/n);
%keyboard

clear Proj

for order=-N:N
  Z=A;
  for i=1:n-1
    P=R^i;
    Z=Z+exp(j*order*2*pi*i/n)*P'*A*P;
  end
  Proj(:,:,order+N+1)=round(1e8*Z/n)/1e8;
  %Proj(:,:,order+N+1)=Z/n;
end

return