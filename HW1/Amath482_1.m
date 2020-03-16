clear; close all; clc;
load Testdata
%1
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1);
x=x2(1:n); 
y=x; 
z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1];
ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
uavg=zeros(64,64,64);
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n);
utemp = fftn(Un);
uavg=uavg+utemp;
end
frq = fftshift(uavg/20);
[mxv,idx] = max(frq(:));
[kxc,kyc,kzc] = ind2sub(size(frq),idx);
KXC = Kx(kxc,kyc,kzc);
KYC = Ky(kxc,kyc,kzc);
KZC = Kz(kxc,kyc,kzc);

%2
tau = 0.3;
filter = exp(-tau*(Kx-KXC).^2).*exp(-tau*(Ky-KYC).^2).*exp(-tau*(Kz-KZC).^2);
filter = fftshift(filter);
trajectory = zeros(20,3);
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n);
Unt = ifftn(fftn(Un).*filter);
[mxv,idx] = max(Unt(:));
[xt,yt,zt] = ind2sub(size(Unt),idx);
trajectory(j,:) = [X(xt,yt,zt) Y(xt,yt,zt) Z(xt,yt,zt)];
end
plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3));
grid on;hold on;
plot3(trajectory(1,1),trajectory(1,2),trajectory(1,3),'r*','MarkerSize',20)
plot3(trajectory(20,1),trajectory(20,2),trajectory(20,3),'.','MarkerSize',20)
xlabel('x')
ylabel('y')
zlabel('z')
title('trajectory of marble')
%3
x20 = trajectory(20,1);
y20 = trajectory(20,2);
z20 = trajectory(20,3);