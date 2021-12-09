clear all;

%% Inputs
mu0 = 4*pi*1e-07;
B0 = mu0;%6e-4;
susc = 1;%0.96; % Magnetic susceptibility
a = 1;%1.4e-6;  % Grain radius, meters
sep=2;    
alpha=0;
L=10;       % Number of multipoles used
debug_mag=1;
%%

f=spherical_harmonic_two_grain(B0,susc, a, sep, alpha, L, debug_mag)
% %%
% sep=2:0.2:4.2;
% susc=0.1;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_0_1(i)=f(3);
% end
% save('perp_sep_dist_0,1.mat','fmag_0_1','sep');
% %%
% sep=2:0.2:4.2;
% susc=0.15;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_0_15(i)=f(3);
% end
% save('perp_sep_dist_0,15.mat','fmag_0_15','sep');
% %%
% sep=2:0.2:4.2;
% susc=0.2;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_0_2(i)=f(3);
% end
% save('perp_sep_dist_0,2.mat','fmag_0_2','sep');
% %%
% sep=2:0.2:4.2;
% susc=0.25;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_0_25(i)=f(3);
% end
% save('perp_sep_dist_0,25.mat','fmag_0_25','sep');
% %%
% sep=2:0.2:4.2;
% susc=0.5;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_0_5(i)=f(3);
% end
% save('perp_sep_dist_0,5.mat','fmag_0_5','sep');
% %%
% sep=2:0.2:4.2;
% susc=0.75;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_0_75(i)=f(3);
% end
% save('perp_sep_dist_0,75.mat','fmag_0_75','sep');
% %%
% sep=2:0.2:4.2;
% susc=1.5;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_1_5(i)=f(3);
% end
% save('perp_sep_dist_1,5.mat','fmag_1_5','sep');
% %%
% sep=2:0.2:4.2;
% susc=2;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_2(i)=f(3);
% end
% save('perp_sep_dist_2.mat','fmag_2','sep');
% %%
% sep=2:0.2:4.2;
% susc=2.5;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_2_5(i)=f(3);
% end
% save('perp_sep_dist_2,5.mat','fmag_2_5','sep');
% %%
% sep=2:0.2:4.2;
% susc=3;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_3(i)=f(3);
% end
% save('perp_sep_dist_3.mat','fmag_3','sep');
% %%
% sep=2:0.2:4.2;
% susc=3.5;
% for i=1:12
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a, sep(i), alpha, L, debug_mag);
%     fmag_3_5(i)=f(3);
% end
% save('perp_sep_dist_3,5.mat','fmag_3_5','sep');
%%
% figure()
% plot(sep,fmag_perp_1(1:12));
% xlabel('sep (a)');
% ylabel('f_z (N)');
% title('H=1 (perpendicular), mu=1, a=1, L=10');
%%
% fprintf('Particle Size Perpedicular 01\n');
% susc=1;
% L=10;
% alpha=90;
% fmag_perp_1=zeros(1,19);
% a=1:0.5:10;
% a=a*1e-6;
% for i=1:19
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a(i), sep, alpha, L);
%     fmag_perp_1(i)=f(3);
% end
% 
% save('particle_size_1_perp.mat','fmag_perp_1','a');
%%
% fprintf('Particle Size Perpedicular 04\n');
% susc=4;
% L=20;
% alpha=90;
% fmag_perp_4=zeros(1,19);
% a=1:0.5:10;
% a=a*1e-6;
% for i=1:19
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0,susc, a(i), sep, alpha, L);
%     fmag_perp_4(i)=f(3);
% end
% 
% save('particle_size_4_perp.mat','fmag_perp_4','a');
%%
% fprintf('B0 variation Parallel 01\n');
% susc=1;
% L=10;
% alpha=0;
% a=1.4e-6;
% fmag_paral_1=zeros(1,19);
% B0=1:0.5:10;
% B0=B0*1e-4;
% for i=1:19
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0(i),susc, a, sep, alpha, L);
%     fmag_paral_1(i)=f(3);
% end
% 
% save('B_mag_1_paral.mat','fmag_paral_1','B0');
%%
% fprintf('B0 variation Parallel 04\n');
% susc=4;
% L=20;
% alpha=0;
% a=1.4e-6;
% fmag_paral_4=zeros(1,19);
% B0=1:0.5:10;
% B0=B0*1e-4;
% for i=1:19
%     fprintf('%i ', i);
%     f=spherical_harmonic_two_grain(B0(i),susc, a, sep, alpha, L);
%     fmag_paral_4(i)=f(3);
% end
% 
% save('B_mag_4_paral.mat','fmag_paral_4','B0');
% %%
% 
% figure(1);
% plot(B0, -fmag_paral_4/mu0);