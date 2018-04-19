%
%  This is a demo of the GLMM algorithm for the Houston Data set.
%
%


clear all;
close all;


rng(5,'twister')

path_var = pwd;
addpath(genpath(path_var))

        
load real_data_1 

v = [57,30,20]; % RGB bands

[m,n,L] = size(data);

rgb(:,:,1) = imadjust(rescale(data(:,:,v(1)),1));
rgb(:,:,2) = imadjust(rescale(data(:,:,v(2)),1));
rgb(:,:,3) = imadjust(rescale(data(:,:,v(3)),1));
figure, imshow(rgb); % display RGB image

data_r = reshape(data,m*n,L);

load endmembers_houston % load initial reference endmembers
R = size(S0,2);
N = m*n;

materials{1} = 'vegetation'; % identify materials
materials{2} = 'red roofs';
materials{3} = 'concrete';
materials{4} = 'asphalt';

r = data_r';

%% Endmember Extraction
% 
% M0 = vca(r,'Endmembers',R);
% id = zeros(R,1);
% for k = 1:R
%     for l = 1:R
%         s(l) = 180*acos( (M(:,k).')*M0(:,l) /(norm(M(:,k))*norm(M0(:,l))) )/pi;
%     end
%     [~, id(k)] = min(s);
% end
% M0 = M0(:,id);

% using endmembers used in [1].
M0 = S0;


%% Fully Constrained Least Squares Unmixing (FCLSU)

disp('FCLSU...')
tic
A_FCLSU = FCLSU(data_r',M0)';
toc
A_FCLSU = A_FCLSU';

R_FCLS = M0*A_FCLSU';

rmse_r_FCLS = RMSEAndSTDForMatrix(r(:), R_FCLS(:));

%% Scaled version of the (partially) Constrained Least Squares (SCLSU)

disp('S-CLSU...')

tic
[A_SCLSU, psis_SCLSU] = SCLSU(data, M0);
toc

S_SCLSU = zeros(L,R,N);

parfor i = 1:m*n
   S_SCLSU(:,:,i) = M0*diag(psis_SCLSU(:,i)); 
end

% S_SCLSU = row2col_lexico_order(S_SCLSU,m,n);
acos_r_SCLS = 0;
acos_r_FCLS = 0;
R_SCLS = zeros(size(r));

for ll=1:N
    R_SCLS(:,ll) = squeeze(S_SCLSU(:,:,ll))*A_SCLSU(:,ll);
    acos_r_SCLS = acos_r_SCLS + acos(r(:,ll)'*R_SCLS(:,ll)/(norm(r(:,ll))*norm(R_SCLS(:,ll))));
    acos_r_FCLS = acos_r_FCLS + acos(r(:,ll)'*R_FCLS(:,ll)/(norm(r(:,ll))*norm(R_FCLS(:,ll))));
end
acos_r_SCLS = acos_r_SCLS/N;
acos_r_FCLS = acos_r_FCLS/N;
rmse_r_SCLS = RMSEAndSTDForMatrix(r(:), R_SCLS(:));




%% Full Generalized Extended Linear Mixing Model

disp('GLMM')


nnorm = '1,1'; % Use a Total Variation on the abundances
verbose = true; % display

maxiter_anls = 50;
maxiter_admm = 100;
epsilon_s = 10^(-3);
epsilon_a = 10^(-3);
epsilon_psi = 10^(-3);
epsilon_admm_abs = 10^(-2);
epsilon_admm_rel = 10^(-2);


% regularization parameters
lambda_m = 1;
lambda_a = 0.01;
lambda_psi = 1e-6;

% initialization of the PSI and Abundance 
psis_init = ones(L,R,N);
A_init = A_SCLSU;

R_GLMM = zeros(size(r));
acos_r_GLMM = 0;


tic
[A_GLMM, psis_GLMM, S_GLMM, optim_struct] = GLMM_ADMM(data, A_init, psis_init, M0,lambda_m,lambda_a,lambda_psi,nnorm,verbose,maxiter_anls,maxiter_admm,epsilon_s,epsilon_a,epsilon_psi,epsilon_admm_abs,epsilon_admm_rel);            
toc

for ll=1:N,
    R_GLMM(:,ll) = squeeze(S_GLMM(:,:,ll))*A_GLMM(:,ll);
    acos_r_GLMM = acos_r_GLMM + acos(r(:,ll)'*R_GLMM(:,ll)/(norm(r(:,ll))*norm(R_GLMM(:,ll))));
end
acos_r_GLMM = acos_r_GLMM/N;
rmse_r_GLMM = RMSEAndSTDForMatrix(r(:), R_GLMM(:));
                        



%%

A_FCLSU_im = reshape(A_FCLSU,m,n,R);
A_SCLSU_im = reshape(A_SCLSU',m,n,R);

A_GLMM_im = reshape(A_GLMM',m,n,R);

fh = figure;
[ha, pos] = tight_subplot(3, 4, 0.01, 0.1, 0.1);
for i=1:4,
    axes(ha(i));
    imagesc(A_FCLSU_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+4));
    imagesc(A_SCLSU_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+8));
    imagesc(A_GLMM_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])
end
% set(fh, 'Position', [0 0 650 700])
axes(ha(1));
ylabel('FCLS','interpreter','latex')
title('Vegetation','interpreter','latex')
axes(ha(2));
title('Met. Roofs','interpreter','latex')
axes(ha(3));
title('Concrete','interpreter','latex')
axes(ha(4));
title('Asphalt','interpreter','latex')
axes(ha(5));
ylabel('SCLSU','interpreter','latex')
axes(ha(9));
ylabel('GLMM','interpreter','latex')

colormap(jet)

fprintf('FCLS: RMSE_R = %f, SAM_R = %f \n', rmse_r_FCLS, acos_r_FCLS)
fprintf('SCLS: RMSE_R = %f, SAM_R = %f \n', rmse_r_SCLS, acos_r_SCLS)
fprintf('GLMM: RMSE_R = %f, SAM_R = %f \n', rmse_r_GLMM,acos_r_GLMM)

