clear all;
%% Inputs

mu0 = 4*pi*1e-07;
B0 = mu0;       %4*pi*1e-07;
susc = 1;       % Magnetic susceptibility
a = 1;          % Grain radius, meters
sep=2;          % Separation between the grains in terms of radius 
alpha=0;        % Magnetic Field Direction
L=10;           % Number of multipoles used
debug_mag=0;    % If 1 plots magnetic field magnitude for each L
debug_f_L=1;    % If 1 plots z component of force with L

mu = (1+susc)*mu0;
H0mag = B0/(mu0);  % Applied magnetic field, A/m
sep=sep*a;

R1=[0 0 sep]';
R2=[0 0 0]';

alpha = deg2rad(alpha); % change angle into radians
H0 =[H0mag*sin(alpha) 0 H0mag*cos(alpha)]'; % A/m

% parallel and perpendicular components of the applied magnetic field
H_perp=H0(1);
H_prll=H0(3);

% %% Solve for coefficients (derived result)
% for m=0:1
%     % vectors for diagonal LXL matrices
%     xm=zeros(1,L);
%     bm=zeros(1,L);
%     cm=zeros(1,L);
%     % q vector
%     qm=zeros(L,1);
%     for l=1:L
%         Xl=(mu/mu0)*l + l + 1;
%         Bl=0;
%         Cl=0;
%         for s=1:L
%             Bl= Bl + ((-1)^(s+m))*nchoosek(l+s,s+m)*(a^(2*l+a)/(sep^(l+s+1)))*((mu/mu0)*l-l);
%             Cl= Cl + ((-1)^(l+s))*((-1)^(s+m))*nchoosek(l+s,s+m)*(a^(2*l+a)/(sep^(l+s+1)))*((mu/mu0)*l-l);
%         end
%         xm(l)=Xl;
%         bm(l)=Bl;
%         cm(l)=Cl;
%     end
%     % LXL matrices
%     Xm=diag(xm);
%     Bm=diag(bm);
%     Cm=diag(cm);
%     %% 2LX2L matrix
%     Am=zeros(2*L);
%     Am(1:L,1:L)=Xm;
%     Am(L+1:end,1:L)=Bm;
%     Am(1:L,L+1:end)=Cm;
%     Am(L+1:end,L+1:end)=Xm;
%     % q vector
%     if m==0
%         qm(1)=-H_prll*(a^3)*(1-mu/mu0);
%     elseif m==1
%         qm(1)=H_perp*(a^3)*(1-mu/mu0);
%     end
%     % 2L Q vector
%     Qm=[qm;qm];
%     %Solve linear equations
%     Beta_m=Am\Qm;
%     Beta1_m=Beta_m(1:L);
%     Beta2_m=Beta_m(L+1:2*L);
%     if m==0
%         Beta1_0=Beta1_m;
%         Beta2_0=Beta2_m;
%     elseif m==1
%         Beta1_1=Beta1_m;
%         Beta2_1=Beta2_m;
%     end
% end

%% Code to solve for coefficients (using paper)
for m=0:1
    % Creating the L X L matrices 
    X=zeros(L,L);
    Delta_m=zeros(L,L);
    Gamma_m=zeros(L,L);
    for i=1:L
        for j=1:L
            % X matrix
            if i==j
                X(i,j)= i*(mu/mu0) + i +1;
            end
            % Delta and Gamma matrix
            Delta_m(i,j)= ((-1)^(i+m))*(i*(mu/mu0)-i)*nchoosek(i+j,j+m)*...
                (a^(2*i+1))/(sep^(i+j+1));
            Gamma_m(i,j)= ((-1)^(i+j))*Delta_m(i,j);
        end
    end
    % 2LX2L matrix
    Am=zeros(2*L);
    Am(1:L,1:L)=X;
    Am(L+1:end,1:L)=Gamma_m;
    Am(1:L,L+1:end)=Delta_m;
    Am(L+1:end,L+1:end)=X;
    
    % q vector
    qm=zeros(L,1);
    if m==0
        qm(1)=-H_prll*(a^3)*(1-mu/mu0);
    elseif m==1
        qm(1)=H_perp*(a^3)*(1-mu/mu0);
    end
    % 2L Q vector
    Qm=[qm;qm];
    %Solve linear equations
    Beta_m=Am\Qm;
    Beta1_m=Beta_m(1:L);
    Beta2_m=Beta_m(L+1:2*L);
    if m==0
        Beta1_0=Beta1_m;
        Beta2_0=Beta2_m;
    elseif m==1
        Beta1_1=Beta1_m;
        Beta2_1=Beta2_m;
    end
end

%% Computing the Magnetic Field

%Create a 3D spherical mesh
dang = pi/18;
inc = dang/2:dang:pi+dang/2;
az = dang/2:dang:2*pi+dang/2;
dr=a/100;
r1=a-dr:dr:a+dr;
[theta,phi,R] = meshgrid(inc,az,r1);

%Convert the spherical mesh to cartesian
x=R.*cos(phi).*sin(theta);
y=R.*sin(phi).*sin(theta);
z=R.*cos(theta);

%define size parameters
size_R=size(R); %size of all elements in 3D space
size_H=size_R; 
size_H(1,4)=L; %size for 4D mag arrays

Hr_L=zeros(size_H);
Hth_L=zeros(size_H);
Hphi_L=zeros(size_H);

Hr=0;
Hth=0;
Hphi=0;

for l=1:L
    for m=0:1
        Hrs=0;
        Hths=0;
        
        for s=m:L
            Psm=legendre(s,cos(theta));
            Ps1m=legendre(s+1,cos(theta));
            if s~=0
                Psm=reshape(Psm(m+1,:,:),size(R));
            end
            Ps1m=reshape(Ps1m(m+1,:,:),size(R));
            %compute derivative of the associated legendre function
            dPsm=((m-s-1).*Ps1m + (s+1).*cos(theta).*Psm)./(-sin(theta));
            % R component
            Hrs=Hrs+ (-1)^(s+m) * nchoosek(l+s,s+m)*s.*(R.^(s-1)).*Psm/...
                (sep^(l+s+1));
            % Theta component
            Hths=Hths+ (-1)^(s+m) * nchoosek(l+s,s+m).*(R.^(s-1)).*dPsm./...
                ((sep^(l+s+1)));
        end
        if m==0
            Plm=legendre(l,cos(theta));
            Pl1m=legendre(l+1,cos(theta));
            Plm=reshape(Plm(m+1,:,:),size(R));
            Pl1m=reshape(Pl1m(m+1,:,:),size(R));
            
            %compute derivative of the associated legendre function
            dPlm=((m-l-1).*Pl1m + (l+1).*cos(theta).*Plm)./(-sin(theta));
            % Theta Component
            Hth=Hth+...
                (Beta1_0(l).*dPlm./((R.^(l+2))) + Beta2_0(l)*Hths).*cos(m*phi);
            % R Component
            Hr=Hr+...
                ((l+1)*Beta1_0(l).*Plm./(R.^(l+2)) -...
                Beta2_0(l)*Hrs).*cos(m*phi);
        elseif m==1
            Plm=legendre(l,cos(theta));
            Pl1m=legendre(l+1,cos(theta));
            Plm=reshape(Plm(m+1,:,:),size(R));
            Pl1m=reshape(Pl1m(m+1,:,:),size(R));
            
            %compute derivative of the associated legendre function
            dPlm=((m-l-1).*Pl1m + (l+1).*cos(theta).*Plm)./(-sin(theta));
            % Theta Component
            Hth=Hth+...
                (Beta1_1(l).*dPlm./((R.^(l+2))) + Beta2_1(l)*Hths).*cos(m*phi);
            % R Component
            Hr=Hr+...
                ((l+1)*Beta1_1(l).*Plm./(R.^(l+2)) -...
                Beta2_1(l)*Hrs).*cos(m*phi);
        end
    end
    Hr_L(:,:,:,l)=Hr;
    Hth_L(:,:,:,l)=-Hth;
end
%Phi component
for l=1:L
    Hphis=0;
    for s=1:L
        Ps1=legendre(s,cos(theta));
        Ps1=reshape(Ps1(2,:,:),size(R));
        Hphis=Hphis+ (-1)^(s+1) * nchoosek(l+s,s+1).*(R.^(s-1)).*Ps1...
            ./(sin(theta).*(sep^(l+s+1)));
    end
    Pl1=legendre(l,cos(theta));
    Pl1=reshape(Pl1(2,:,:),size(R));
    Hphi=Hphi+...
                (Beta1_1(l).*Pl1./(sin(theta).*(R.^(l+2))) +...
                Beta2_1(l)*Hphis).*sin(phi);

    Hphi_L(:,:,:,l)=Hphi;
end
Hth=-Hth;

%% plot Magnetic Field
if debug_mag==1
    for l=1:L
        H_tot_mag=zeros(size_R);
        for i=1:size_R(1)
            for j=1:size_R(2)
                for k=1:size_R(3)
                    ph=az(i);
                    th=inc(j);
                    %H_part_sph=[Hr(i,j,k);Hphi(i,j,k);Hth(i,j,k)];
                    %transformation matrix
                    pre=[sin(th)*cos(ph), cos(th)*cos(ph), -sin(ph);...
                        sin(th)*sin(ph), cos(th)*sin(ph),  cos(ph);...
                        cos(th),         -sin(th),         0];
                    post=pre';
                    H0_sph=post*H0;
                    H_sph=[Hr_L(i,j,k,l);Hth_L(i,j,k,l);Hphi_L(i,j,k,l)] + H0_sph;
                    H_tot_mag(i,j,k)=norm(H_sph);
                end
            end
        end
        figure();
        Ang=1;
        pc=pcolor(squeeze(z(Ang,:,:))./a,squeeze(x(Ang,:,:))./a,squeeze(H_tot_mag(Ang,:,:))); set(pc, 'EdgeColor', 'none');
        colormap('hot');
        xlim([-2 2]); ylim([-1.5 1.5]);
        title('|H|');
        xlabel('x');
        ylabel('y');
        colorbar;
        axis equal;
        grid on;
    end
end

%% Formulating the Maxwell Stress Tensor in Spherical Coordinates
f_L=zeros(3,1,L);
for l=1:L
    f=0;
    for i=1:size_R(1)
        if i==1 || i==size_R(1)
            p=1;
        elseif mod(i,2)~=0
            p=2;
        else
            p=4;
        end
        for j=1:size_R(2)
            if j==1 || j==size_R(2)
                q=1;
            elseif mod(j,2)~=0
                q=2;
            else
                q=4;
            end
            ph=az(i);
            th=inc(j);
            %transformation matrix
            pre=[sin(th)*cos(ph), cos(th)*cos(ph), -sin(ph);...
                sin(th)*sin(ph), cos(th)*sin(ph),  cos(ph);...
                cos(th),         -sin(th),         0];
            post=pre';
            H0_sph=post*H0;
            H_sph=[Hr_L(i,j,2,l);Hth_L(i,j,2,l);Hphi_L(i,j,2,l)] + H0_sph;
            H_cart=pre*H_sph;
            h=norm(H_cart);
            T_cart=mu0*(H_cart*H_cart' - 0.5*h^2*eye(3,3));
            rn_hat=[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]';
            f=f+a*a*T_cart*rn_hat*sin(th)*p*q;
        end          
    end
    f=f*dang*dang/9;
    f_L(:,:,l)=f;
end

if debug_f_L==1
    figure()
    plot(1:L,squeeze(f_L(3,1,:)))
    title('Z component of force vs L (multipoles)')
    xlabel('L')
    ylabel('f_z (N)')
    grid on
end