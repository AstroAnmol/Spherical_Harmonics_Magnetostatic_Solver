clear all;

%% Inputs
mu0 = 4*pi*1e-07;
B0 = 6e-4;
susc = 0.96; % Magnetic susceptibility
a = 1.4e-6;  % Grain radius, meters
sep=2.2;
alpha=90;
L=10;       %Number of multipoles used

%%

f=spherical_harmonic_two_grain(B0,susc, a, sep, alpha, L)
%%
fmag=zeros(1,12);
sep=2:0.2:4.2;
for i=1:12
    fprintf('%i ',i);
    f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L);
    fmag(i)=f(3);
end

save('keaveny_1_parallel.mat');
%%
figure(1);
plot(sep, -fmag/mu0);


% l=5;
% theta=pi/4;
% m=3;
% Plm1=legendre(l-1,cos(theta));
% Plm=legendre(l,cos(theta));
% Pl1m=legendre(l+1,cos(theta));
% Plm=reshape(Plm(m+1,:,:),size(theta));
% Pl1m=reshape(Pl1m(m+1,:,:),size(theta));
% Plm1=reshape(Plm1(m+1,:,:),size(theta));
%             
% %compute derivative of the associated legendre function
% dPlm=((m-l-1).*Pl1m + (l+1).*cos(theta).*Plm)./(-sin(theta));
% 
% DPlm=(l*cos(theta).*Plm - (l+m)*Plm1)/(sqrt(1-cos(theta).*cos(theta)));