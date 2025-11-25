%% Script description

% This demo script can reconstruct data acquired from the 4-shot mosaic sequences,
% in ...\Pulseq_sequence\,
% tested with a Siemens 3T Prisma scanner at the Max Planck Institute for Biological Cybernetics in Tuebingen, Germany.

% For questions, please contact: Rui Tian, rui.tian@tuebingen.mpg.de/ruitianwater@outlook.com

% Please cite the following conference/workshop abstracts, and the incoming research articles for this work:
% Tian R, Uecker M, Zaitsev M, Scheffler K. Overlap-kernel EPI: estimating MRI shot-to-shot phase variations by shifted-kernel extraction from overlapped regions at arbitrary k-space locations. Magn Reson Med, 2025. doi: 10.1002/mrm.70196
% Tian R, Uecker M, Scheffler K. Navigator-free multi-shot EPI with shift-invariant kernel extraction in subspace. Power pitch. Program #1081. In proceedings of the annual meeting of ISMRM, Singapore, 2024.
% Tian R, Uecker M, Zaitsev M, Scheffler K. Pulseq Implementation of Overlapped Readout Segmented EPI with Phase Fluctuations Corrected by Shift-Invariant Kernel Extraction. Traditional poster #18. ISMRM Workshop on 40 Years of Diffusion: Past, Present & Future Perspectives, Kyoto, Japan, 2025.

% add subfolders
addpath(genpath('./'))

% MATLAB stops at error locations
dbstop error;

tst_tot = tic;

%% Readme

% The current protocol:
% 4-shot, mosaic segmented, 2-slice, b=0,1000s/mm2, 2D navigator ON

% Notice: This sequence is more sensitive to B0 inhomogeneity, compared to
% readout-segmented and phase-interleaved schemes. Thus, the 1st and the
% 2nd slices of this sample data show reconstructions with and without
% noticable distortion artifacts, respectively. Further improvements of
% distortion corrections are expected.

% Additionally, this 4-shot sequence with self-navigated inter-shot correction
% should apparently reduce shot-to-shot eddy currents. However, residue
% artifacts due to trajectory errors/even odd echoes corrections might be
% observed but beyond the invesigation of this paper. 
% This should be investigated in future work.

%% User input: enter sequence/recon parameters

% Experiment number
para.nEPI = 224; % EPI scan MID
para.nNS = 218; % Noise scan MID

% Intershot correction for multi-shot EPI
para.bIntSCorr=1; % 1: Apply inter-segment correction, 0: does not apply
para.type_subspace = 1; % 1: ESPIRiT 2: GRAPPA
para.nky_subspace=7; % kernel length along y
para.nkx_subspace=7; % kernel length along x
para.thres_k_diff = round(para.nky_subspace*para.nkx_subspace*2*0.5); % empirical kernel threshold, b=0 s/mm2
para.thres_k_b0 = round(para.nky_subspace*para.nkx_subspace*2*0.5); % empirical kernel threshold, b~=0 s/mm2

% Parallel imaging
para.bPwhite = 1; % Pre-whitening, if noise data is available
para.iterLSQR_sepShot = 5; % Max. iteration for LSQR, intermediate PI recon.
para.iterLSQR_combShot = 10; % Max. iteration for LSQR, combining shots

% 2D Navigator
para.bNav = 0; % 1: Use 2D navigator 0: disable 2D navigator. 
para.bConjNav =0; % take conjugate of navigator images before applying, 1 for GE-navigator

% Common sequence parameter
para.bAmpCorr = 0; % 1: shot-to-shot amplitude correction included. 0: otherwise
para.nOS = 1; % readout oversampling factor
para.nPI = 3; % parallel imaging acceleration factor
para.nPFT = 8/8; % partial Fourier
para.coeff_kai = 3; % Kaiser filter coeff, short readout axis smoothing
para.bMolComb = 0; % 0: sum up shot images 1: combine shots in a big forward model, be careful about the convergence.
para.bDisCorr = 1; % 1: Use distortion correction, 0: no distortion correction.
para.deltaTE = 3; % For 2 GRE field mapping scans, TE2-TE1 in ms
para.MID_FieldMap = [220,221]; % The MIDs for GRE scans with 2 echo times.
para.nCali = 219; % GRE calibration scan for parallel imaging
para.kB0 = 20; % kernel window size of B0 mapping. Try smaller values (e.g., 10) for reducing recon instability given too strong distortion.

% Scan dependent sequence parameter
switch para.nEPI
    case 224
        para.seq_file_path = 'rtian_03_4mo_SE_epi_b1k.seq'; % multi-shot EPI sequence file
        nProj = 229; % MID, 1D navigator scan for odd-even echo correction
        para.seq_proj_file_path = 'rtian_03_proj.seq';  % sequence file for 1D navigator scan
        para.ESP = 0.7; % ESP [ms], see printout of the Pulseq sequence
        para.ind_kc = 4; % index of the central k-space phase encoded lines, for distortion corr. in SE-EPI
        para.ind_Nav_kc = 7; % index of the central k-space phase encoded lines, for distortion corr. in SE-EPI
end

if (para.bNav==1) && (para.bIntSCorr ==1)
    error('Please only choose one inter-shot correction method.')
end

%% Below: No user input needed

%% Load Siemens data

% Load multishot EPI data
file_epi_ms = mapVBVD(para.nEPI);
file_epi_ms = getLastStruct(file_epi_ms);
dat_epi = file_epi_ms.image();
dat_epi = permute(dat_epi,[1,2,3,4,10,6,7,8,9,5,11]); % Switch SET and SLC label, used to avoid pulseq slice limitations

p = size(dat_epi,5);
new_order = reshape([1:2:p,2:2:p],1,[]);
perm_order(new_order) = 1:p;
dat_epi = dat_epi(:,:,:,:,perm_order,:,:,:,:,:,:);

% Load navigator data, if needed
if para.bNav==1
    dat_2DNav = file_epi_ms.phasecor();
    dat_2DNav = permute(dat_2DNav,[1,2,3,4,10,6,7,8,9,5,11]);% Switch SET and SLC label, used to avoid pulseq s
    dat_2DNav = dat_2DNav(:,:,:,:,perm_order,:,:,:,:,:,:);
    dat_2DNav(:,:,:,:,:,:,:,:,:,:,[1,3]) = flip(dat_2DNav(:,:,:,:,:,:,:,:,:,:,[1,3]),3); % for 4s-EPI, the 2D navigators have up-down directions flip
end

% Size of EPI data [Read,RF,Phase,1,2D Slice,1,1,1,1,1,Segment/Shot]
[para.NCol,para.NCha,para.NLin,para.NPar,para.NSli,para.NAve,para.NPhs,para.NEco,para.NRep,para.NSet,para.NSeg] = size(dat_epi);

%% Load noise scan

if para.bPwhite ==1
    noise = mapVBVD(para.nNS).noise();
    noise = permute(noise,[2,1,3,4,5,6]);
    noise = reshape(noise,size(noise,1),[]);
    mat_noise = (noise*noise')./(size(noise,2)-1);
    L = chol(mat_noise,'lower');
    L_inv = inv(L);
end

%% Re-grid EPI acquisition

% Load trajectory from .seq files
seq = mr.Sequence(); % Create a new sequence object
seq.read(para.seq_file_path,'detectRFuse');
[ktraj_adc, t_adc, ktraj, ~, t_excitation, t_refocusing] = seq.calculateKspacePP();
k_acq = ktraj_adc;
k_acq = reshape(k_acq,3,para.NCol,[],para.NPar,para.NSli,para.NSeg,para.NPhs,para.NRep,para.NAve,para.NEco,para.NSet);
k_acq = permute(k_acq,[1,2,3,4,5,9,7,10,8,11,6]);

% Extract 2D navigator trajectory
if para.bNav==1
    para.NLin_2DNav = size(k_acq,3)-para.NLin;
    k_2DNav = k_acq(:,:,end-para.NLin_2DNav+1:end,:,:,:,:,:,:,:,:);
end

% imaging trajectory
k_acq = k_acq(:,:,1:para.NLin,:,:,:,:,:,:,:,:);

% EPI gridding
[dat_epi_rg,P1,P2,para.mapS] = grid2DEPI(dat_epi,k_acq,para.nOS,para.nPI,1,para);
para.ny_tot = size(dat_epi_rg,1);
para.nx_tot = length(P1);

% 2D navigator gridding
if para.bNav==1
    [dat_2DNav_rg,P_Nav,para.mapS_Nav] = gridEPI(dat_2DNav,k_2DNav,para.nOS,para.nPI,2,para);
    ny_Nav = size(dat_2DNav_rg,1);
    nx_Nav = length(P_Nav);
end

% k-space coverage
ref_MTF_shots = double(squeeze((dat_epi_rg(:,1,1,1,1,1,1,1,1,1,:)~=0)));
ref_MTF = sum(ref_MTF_shots,2)./2;
ref_PE_MTF = (P1==1)&(P2==1);
ind_MTF = find(ref_PE_MTF==1);
ref_PE_MTF((ind_MTF(1)):(ind_MTF(end)))=1;

% k-space masks for 4 bands. The codes can be optimized in future work.
ref_MTF_PE_temp = zeros(1,length(P1));
ref_MTF_PE_temp(1:ind_MTF(end)) = 1;
M1 = ref_MTF_shots(:,1)*ref_MTF_PE_temp;
K1 = kaiser(sum(ref_MTF_shots(:,1)==1),para.coeff_kai)...
    *kaiser(sum(ref_MTF_PE_temp==1),para.coeff_kai)';
K1_temp= M1;K1_temp(K1_temp==1)=K1;
K1 = K1_temp;

ref_MTF_PE_temp = zeros(1,length(P1));
ref_MTF_PE_temp(ind_MTF(1):end) = 1;
M2 = ref_MTF_shots(:,2)*ref_MTF_PE_temp;
K2 = kaiser(sum(ref_MTF_shots(:,2)==1),para.coeff_kai)...
    *kaiser(sum(ref_MTF_PE_temp==1),para.coeff_kai)';
K2_temp= M2;K2_temp(K2_temp==1)=K2;
K2 = K2_temp;

ref_MTF_PE_temp = zeros(1,length(P1));
ref_MTF_PE_temp(1:ind_MTF(end)) = 1;
M3 = ref_MTF_shots(:,3)*ref_MTF_PE_temp;
K3 = kaiser(sum(ref_MTF_shots(:,3)==1),para.coeff_kai)...
    *kaiser(sum(ref_MTF_PE_temp==1),para.coeff_kai)';
K3_temp= M3;K3_temp(K3_temp==1)=K3;
K3 = K3_temp;

ref_MTF_PE_temp = zeros(1,length(P1));
ref_MTF_PE_temp(ind_MTF(1):end) = 1;
M4 = ref_MTF_shots(:,4)*ref_MTF_PE_temp;
K4 = kaiser(sum(ref_MTF_shots(:,4)==1),para.coeff_kai)...
    *kaiser(sum(ref_MTF_PE_temp==1),para.coeff_kai)';
K4_temp= M4;K4_temp(K4_temp==1)=K4;
K4 = K4_temp;

M = cat(3,M1,M2,M3,M4);
K_2D = cat(3,K1,K2,K3,K4);
M_sum = sum(M,3);

dat_epi = dat_epi_rg; clear dat_epi_rg; % clear to reduce memory
dat_epi = permute(dat_epi,[1,3,11,2,5,7,4,6,8,9,10]); %-> update dim. [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]

% data for 2D navigator
if para.bNav==1
    dat_2DNav = permute(dat_2DNav_rg,[1,3,11,2,5,7,4,6,8,9,10]); %-> update dim. [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    clear dat_2DNav_rg;
end

disp('Naive EPI acquisition gridding done')

%% Load 1D navigator scan, odd/even echoes for N/2 ghost correction

% Load 1D projection line
file_proj = mapVBVD(nProj);
file_proj = getLastStruct(file_proj);
dat_proj = file_proj.image();

% The size of projection line
sz_proj = size(dat_proj);

%% Regrid odd/even echoes

% Load trajectory from .seq files
seq = mr.Sequence(); % Create a new sequence object
seq.read(para.seq_proj_file_path,'detectRFuse');
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
k_proj = ktraj_adc;
k_proj = reshape(k_proj,3,sz_proj(1),[]);

[dat_proj_rg,~,~] = gridEPI(dat_proj,k_proj,para.nOS,0,1,para);

dat_proj=dat_proj_rg; clear dat_proj_rg; % clear to reduce memory
ny_echoes = size(dat_proj,1);

figure
for idx_fig = 1:sz_proj(2)
    subplot(round(sqrt(sz_proj(2))),round(sqrt(sz_proj(2))),idx_fig)
    imagesc(abs(squeeze(dat_proj(:,idx_fig,:))))
    title(['projection line, coil#',num2str(idx_fig)])
end

disp('Naive regridding of projection data done.')

%% Correct EPI ghost with the 1D projection line

% sum up k-space data from RF channels
line_proj = SymFft(dat_proj,1);

% Steps of phase encoded lines to take into account, empirical
nSteps = 15;

% Calculate difference between neighboring lines
offset = zeros(size(line_proj,1),size(line_proj,2),2,length(2:2:(2*nSteps)));

ct = 0;
for idx_PE = 2:2:(2*nSteps)

    ct = ct+1;

    p1 = line_proj(:,:,idx_PE-1,:)./...
        line_proj(:,:,idx_PE,:);
    p2 = line_proj(:,:,idx_PE+1,:)./...
        line_proj(:,:,idx_PE,:);

    offset(:,:,1,ct) = p1;
    offset(:,:,2,ct) = p2;
end

% Average the difference, including RF receivers dim.
offset_avg = mean(real(offset),[2 3 4]) + 1i.*mean(imag(offset),[2 3 4]);

% Plot the result
figure
plot(angle(offset_avg(:,:,1,1,1)))
title('phase ramp, even odd echoes')

% Apply the N/2 ghost correction
bIntpHR = 1; % interpolate the even/odd echoes shift to high resolution
b4shot = 1;%rtian-4shot
[dat_epi_corr,img_epi_1Dcorr] = corr_EvenOddEchoes(dat_epi,offset_avg,ref_MTF_shots,para.ny_tot,ny_echoes,bIntpHR,3*(b4shot==1)); %rtian-4shot
if para.bNav==1
    [dat_nav_1Dcorr,img_nav_1Dcorr] = corr_EvenOddEchoes(dat_2DNav,offset_avg,ref_MTF_shots,ny_echoes,ny_echoes,bIntpHR,0);
end

% the N/2 ghost corrected image
img_epi_1Dcorr_sos = squeeze(sqrt(sum(img_epi_1Dcorr.*conj(img_epi_1Dcorr),4)));

disp('Naive N/2 ghost correction done.')

%% Pre-process data

% Transform the k-space dat to image space
img_epi = SymFft(dat_epi_corr,[1 2]);

% start from multi segment k-space data
dat_mb = dat_epi_corr;
loc_pe = find(P1==1);loc_pe=loc_pe(1); % for partial Fourier only

% 2D navigator
if para.bNav==1
    dat_nav_mb = dat_nav_1Dcorr;
end

%% Reconstruct EPI images with SENSE

if para.nPI>=1

    % margin recon
    n_mar = 0;

    % ESPIRiT parameters
    nk = 7;
    ny = para.ny_tot+n_mar*2;
    nx = para.nx_tot+n_mar*2;
    thres_k = 0.0121;

    % Load the GRE images for estimating SENSE maps
    file_cali = mapVBVD(para.nCali);
    file_cali = getLastStruct(file_cali);
    dat_cali = file_cali.image.unsorted();
    sz_cali = size(dat_cali);
    dat_cali = reshape(dat_cali,sz_cali(1),sz_cali(2),[],para.NSli);

    % The GRE sos image
    img_cali = SymFft(dat_cali,[1 3]);
    img_cali_sos = squeeze(sqrt(sum(img_cali.*conj(img_cali),2)));

    % ESPIRiT RF maps estimation
    dat_cali = permute(dat_cali,[1 3 2 4]);
    if exist(para.mapS,'file') ~= 2
        maps_sense = zeros(ny,nx,size(dat_cali,3),para.NSli);
        maps_eigval = zeros(ny,nx,size(dat_cali,3),para.NSli);
        for idx_5 = 1:para.NSli
            lowRes_dat = dat_cali(:,:,:,idx_5);
            [eigen_maps,eigen_vals] = ESPIRiT_map(lowRes_dat,nk,ny,nx,thres_k);
            maps_sense(:,:,:,idx_5) = eigen_maps(:,:,:,1);
            maps_eigval(:,:,:,idx_5) = eigen_vals(:,:,:,1);
        end
        maps_mask = maps_eigval(:,:,1,:)>(0.99.*max(maps_eigval(:,:,1,:),[],[1,2]));
        save(para.mapS,'maps_sense','maps_eigval','maps_mask')
    else
        load(para.mapS);
    end

    % rtian-4shot: Load the GRE images for estimating B0 maps
    if para.bDisCorr==1
        para.mapB0 = ['maps_B0_',num2str(para.ny_tot),'_',num2str(para.nx_tot),'_k',num2str(para.kB0),'_MID',num2str(para.nEPI),'.mat'];
        if exist(para.mapB0,'file') ==2
            load(para.mapB0);
        else
            file_B0map_TE1 = mapVBVD(para.MID_FieldMap(1));
            file_B0map_TE1 = getLastStruct(file_B0map_TE1);
            dat_B0map_TE1 = file_B0map_TE1.image.unsorted();

            sz_B0map_TE1 = size(dat_B0map_TE1);
            dat_B0map_TE1 = reshape(dat_B0map_TE1,sz_B0map_TE1(1),sz_B0map_TE1(2),[],para.NSli);

            file_B0map_TE2 = mapVBVD(para.MID_FieldMap(2));
            file_B0map_TE2 = getLastStruct(file_B0map_TE2);
            dat_B0map_TE2 = file_B0map_TE2.image.unsorted();

            sz_B0map_TE2 = size(dat_B0map_TE2);
            dat_B0map_TE2 = reshape(dat_B0map_TE2,sz_B0map_TE1(1),sz_B0map_TE1(2),[],para.NSli);

            B0_fmaps = zeros(ny,nx,2,para.NSli);
            for idx_sl = 1:para.NSli
                dat_B0maps_2TE = cat(4,dat_B0map_TE1(:,:,:,idx_sl),dat_B0map_TE2(:,:,:,idx_sl));
                dat_B0maps_2TE = reshape(dat_B0maps_2TE,[],size(dat_B0maps_2TE,3),1,size(dat_B0maps_2TE,4));
                [B0_Emaps,~] = ESPIRiT_B0map(dat_B0maps_2TE,para.kB0,para.kB0,ny,nx,para.kB0*para.kB0*2*0.5);
                B0_fmaps(:,:,:,idx_sl) = B0_Emaps(:,:,:,1);
            end
            B0_fmaps = squeeze(B0_fmaps(:,:,2,:)./B0_fmaps(:,:,1,:));
            save(para.mapB0,'B0_fmaps');
        end
    end

    % ESPIRiT RF maps estimation for 2D navigator's size, can also be resized
    % from previous estimations.
    if para.bNav==1
        if exist(para.mapS_Nav,'file') ~=2
            maps_Nav_sense = zeros(ny_Nav,nx_Nav,size(dat_cali,3),para.NSli);
            maps_Nav_eigval = zeros(ny_Nav,nx_Nav,size(dat_cali,3),para.NSli);
            parfor idx_5 = 1:para.NSli
                lowRes_dat = dat_cali(:,:,:,idx_5);
                [eigen_maps,eigen_vals] = ESPIRiT_map(lowRes_dat,nk,ny_Nav,nx_Nav,thres_k);
                maps_Nav_sense(:,:,:,idx_5) = eigen_maps(:,:,:,1);
                maps_Nav_eigval(:,:,:,idx_5) = eigen_vals(:,:,:,1);
            end
            maps_Nav_mask = maps_Nav_eigval(:,:,1,:)>(0.9990.*max(maps_Nav_eigval(:,:,1,:),[],[1,2]));
            save(para.mapS_Nav,'maps_Nav_sense','maps_Nav_eigval','maps_Nav_mask')
        else
            load(para.mapS_Nav);
        end
    end

    disp('ESPIRiT maps done')

    % pre-whitening
    if para.bPwhite==1
        dat_mb = permute(dat_mb,[4,1,2,3,5,6,7,8,9,10,11]);
        dat_mb = pagemtimes(L_inv,dat_mb);
        dat_mb = permute(dat_mb,[2,3,4,1,5,6,7,8,9,10,11]);
        maps_sense = permute(maps_sense,[3,1,2,4]);
        maps_sense = pagemtimes(L_inv,maps_sense);
        maps_sense = permute(maps_sense,[2,3,1,4]);
        if para.bNav ==1
            dat_nav_mb = permute(dat_nav_mb,[4,1,2,3,5,6,7,8,9,10,11]);
            dat_nav_mb = pagemtimes(L_inv,dat_nav_mb);
            dat_nav_mb = permute(dat_nav_mb,[2,3,4,1,5,6,7,8,9,10,11]);
            maps_Nav_sense = permute(maps_Nav_sense,[3,1,2,4]);
            maps_Nav_sense = pagemtimes(L_inv,maps_Nav_sense);
            maps_Nav_sense = permute(maps_Nav_sense,[2,3,1,4]);
        end
    end

    % Below: Apply SENSE reconstruction
    dat_mb = padarray(dat_mb,[n_mar 0 0 0 0 0 0 0 0 0 0],0,'both');
    dat_mb=dat_mb.*reshape(ref_MTF_shots,[],1,size(ref_MTF_shots,2));
    dat_mb_y_kx = SymFft(dat_mb,1);
    if para.bNav==1
        dat_nav_mb_y_kx = SymFft(dat_nav_mb,1);
    end

    % Construct Fourier matrix
    E = dftmtx(para.nx_tot+n_mar*2);
    E = circshift(conj(E),[floor(para.nx_tot/2+n_mar) floor(para.nx_tot/2+n_mar)]);
    E = E((n_mar+1):(end-n_mar),:);
    if para.bNav==1
        E_Nav = dftmtx(nx_Nav);
        E_Nav = circshift(conj(E_Nav),[floor(nx_Nav/2) floor(nx_Nav/2)]);
    end

    % Undersampled k-space
    E_P1 = E(P1==1,:);
    E_P2 = E(P2==1,:);E_P2 = flip(E_P2,1);

    % rtian-2shot, empirical and works for scanner data dimensions
    E_P1_temp = E_P1;
    E_P1 = E_P2;
    E_P2 = E_P1_temp;
    clear E_P1_temp;

    % Expand the encoding matrix for sense
    E_P1 = repmat(E_P1,[para.NCha,1]);
    E_P2 = repmat(E_P2,[para.NCha,1]);
    if para.bNav==1
        E_Nav = E_Nav(P_Nav==1,:);
        E_Nav = repmat(E_Nav,[para.NCha,1]);
    end

    % dat_mb: [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    % img_sense: [Read,Phase,Seg,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    img_sense = zeros(ny,nx,size(dat_mb_y_kx,3),size(dat_mb_y_kx,5),...
        size(dat_mb_y_kx,6),size(dat_mb_y_kx,7),size(dat_mb_y_kx,8),...
        size(dat_mb_y_kx,9),size(dat_mb_y_kx,10),size(dat_mb_y_kx,11));

    if para.bNav==1
        img_Nav = zeros(ny_Nav,nx_Nav,1,size(dat_nav_mb_y_kx,5),...
            size(dat_nav_mb_y_kx,6),size(dat_nav_mb_y_kx,7),size(dat_nav_mb_y_kx,8),...
            size(dat_nav_mb_y_kx,9),size(dat_nav_mb_y_kx,10),size(dat_nav_mb_y_kx,11));
    end

    n3_temp  = size(dat_mb_y_kx,3);
    n5_temp  = size(dat_mb_y_kx,5);
    n6_temp  = size(dat_mb_y_kx,6);
    n7_temp  = size(dat_mb_y_kx,7);
    n8_temp  = size(dat_mb_y_kx,8);
    n9_temp  = size(dat_mb_y_kx,9);
    n10_temp  = size(dat_mb_y_kx,10);
    n11_temp  = size(dat_mb_y_kx,11);

    if para.bDisCorr==1
        B0_fmaps_temp = B0_fmaps;
    else
        B0_fmaps_temp=1;
    end

    parfor idx_RO = 1:ny
        for idx_3 = 1:n3_temp
            for idx_5 = 1:n5_temp
                for idx_6 = 1:n6_temp
                    for idx_7 = 1:n7_temp
                        for idx_8 = 1:n8_temp
                            for idx_9 = 1:n9_temp
                                for idx_10 = 1:n10_temp
                                    for idx_11 = 1:n11_temp

                                        dat_temp = dat_mb_y_kx(idx_RO,:,idx_3,:,idx_5,...
                                            idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);
                                        dat_temp = dat_temp(:);

                                        temp_sense = maps_sense(idx_RO,:,:,idx_5);
                                        temp_sense = repmat(temp_sense,[sum(P1==1),1,1]);
                                        temp_sense = permute(temp_sense,[1,3,2]);
                                        temp_sense = reshape(temp_sense,[],nx);

                                        if para.bDisCorr==1
                                            B0maps_temp = B0_fmaps_temp(idx_RO,:,idx_5);
                                            B0_scale = (-para.ind_kc+(1:sum(P1==1)))./para.deltaTE.*para.ESP;B0_scale = reshape(B0_scale,[],1);
                                            B0maps_temp = exp(1i.*angle(B0maps_temp).*B0_scale);
                                            B0maps_temp = repmat(B0maps_temp,[para.NCha,1]);
                                            B0maps_temp = B0maps_temp;
                                        else
                                            B0maps_temp=B0_fmaps_temp;
                                        end

                                        if (idx_3==1) || (idx_3==3)
                                            E_sense = E_P1.*temp_sense.*B0maps_temp;
                                        elseif (idx_3==2) || (idx_3==4)
                                            E_sense = E_P2.*temp_sense.*B0maps_temp;
                                        else
                                            error('error: only 2 bands are allowed')
                                        end

                                        [x_temp,flag,relres,iter,resvec,lsvec] = lsqr(E_sense,dat_temp,[],para.iterLSQR_sepShot);
                                        
                                        % Or other solvers, like pseudoinverse
                                        % x_temp = pinv(E_sense)*dat_temp;
                                        
                                        img_sense(idx_RO,:,idx_3,idx_5,...
                                            idx_6,idx_7,idx_8,idx_9,idx_10,idx_11) = x_temp;

                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    clear dat_mb_y_kx; clear dat_temp; % clear to reduce memory
    img_sense = SymIfft(img_sense,[1 2]);
    img_sense = img_sense.*M;
    img_sense = SymFft(img_sense,[1 2]);

    if para.bNav==1
        % B0 maps for 2D navigator'size, can also be resized from previous
        % estimations.
        if para.bDisCorr==1
            para.mapB0_Nav = ['maps_Nav_B0_',num2str(para.ny_tot),'_',num2str(para.nx_tot),'_k',num2str(para.kB0),'_MID',num2str(para.nEPI),'.mat'];
            if exist(para.mapB0_Nav,'file') ==2
                load(para.mapB0_Nav);
            else
                file_B0map_TE1 = mapVBVD(para.MID_FieldMap(1));
                file_B0map_TE1 = getLastStruct(file_B0map_TE1);
                dat_B0map_TE1 = file_B0map_TE1.image.unsorted();

                sz_B0map_TE1 = size(dat_B0map_TE1);
                dat_B0map_TE1 = reshape(dat_B0map_TE1,sz_B0map_TE1(1),sz_B0map_TE1(2),[],para.NSli);

                file_B0map_TE2 = mapVBVD(para.MID_FieldMap(2));
                file_B0map_TE2 = getLastStruct(file_B0map_TE2);
                dat_B0map_TE2 = file_B0map_TE2.image.unsorted();

                sz_B0map_TE2 = size(dat_B0map_TE2);
                dat_B0map_TE2 = reshape(dat_B0map_TE2,sz_B0map_TE1(1),sz_B0map_TE1(2),[],para.NSli);

                B0_Nav_fmaps = zeros(para.ny_tot,length(P_Nav),2,para.NSli);
                for idx_sl = 1:para.NSli
                    dat_B0maps_2TE = cat(4,dat_B0map_TE1(:,:,:,idx_sl),dat_B0map_TE2(:,:,:,idx_sl));
                    dat_B0maps_2TE = reshape(dat_B0maps_2TE,[],size(dat_B0maps_2TE,3),1,size(dat_B0maps_2TE,4));
                    [B0_Emaps,~] = ESPIRiT_B0map(dat_B0maps_2TE,para.kB0,para.kB0,para.ny_tot,length(P_Nav),para.kB0*para.kB0*2*0.5);
                    B0_Nav_fmaps(:,:,:,idx_sl) = B0_Emaps(:,:,:,1);
                end
                B0_Nav_fmaps = squeeze(B0_Nav_fmaps(:,:,2,:)./B0_Nav_fmaps(:,:,1,:));
                save(para.mapB0_Nav,'B0_Nav_fmaps');
            end
        end

        n3_temp  = size(dat_nav_mb_y_kx,3);
        n5_temp  = size(dat_nav_mb_y_kx,5);
        n6_temp  = size(dat_nav_mb_y_kx,6);
        n7_temp  = size(dat_nav_mb_y_kx,7);
        n8_temp  = size(dat_nav_mb_y_kx,8);
        n9_temp  = size(dat_nav_mb_y_kx,9);
        n10_temp  = size(dat_nav_mb_y_kx,10);
        n11_temp  = size(dat_nav_mb_y_kx,11);

        if para.bDisCorr==1
            B0_Nav_fmaps_temp = B0_Nav_fmaps;
        else
            B0_Nav_fmaps_temp=1;
        end

        % SENSE reconstruction of 2D navigator data
        parfor idx_RO = 1:ny_Nav
            for idx_3 = 1:n3_temp
                for idx_5 = 1:n5_temp
                    for idx_6 = 1:n6_temp
                        for idx_7 = 1:n7_temp
                            for idx_8 = 1:n8_temp
                                for idx_9 = 1:n9_temp
                                    for idx_10 = 1:n10_temp
                                        for idx_11 = 1:n11_temp

                                            dat_nav_temp = dat_nav_mb_y_kx(idx_RO,:,idx_3,:,idx_5,...
                                                idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);
                                            dat_nav_temp = dat_nav_temp(:);

                                            temp_sense = maps_Nav_sense(idx_RO,:,:,idx_5);
                                            temp_sense = repmat(temp_sense,[sum(P_Nav==1),1,1]);
                                            temp_sense = permute(temp_sense,[1,3,2]);
                                            temp_sense = reshape(temp_sense,[],nx_Nav);

                                            if para.bDisCorr==1
                                                B0maps_temp = B0_Nav_fmaps_temp(idx_RO,:,idx_5);
                                                B0_scale = (-para.ind_Nav_kc+(1:sum(P_Nav==1)))./para.deltaTE.*para.ESP;B0_scale = reshape(B0_scale,[],1);
                                                B0maps_temp = exp(1i.*angle(B0maps_temp).*B0_scale);
                                                B0maps_temp = repmat(B0maps_temp,[para.NCha,1]);
                                            else
                                                B0maps_temp=B0_Nav_fmaps_temp;
                                            end


                                            E_sense = E_Nav.*temp_sense.*B0maps_temp;

                                            [x_temp,flag,relres,iter,resvec,lsvec] = lsqr(E_sense,dat_nav_temp);

                                            % Or, other solvers like
                                            % pseudoinverse
                                            % x_temp = pinv(E_sense)*dat_temp;

                                            img_Nav(idx_RO,:,idx_3,idx_5,...
                                                idx_6,idx_7,idx_8,idx_9,idx_10,idx_11) = x_temp;

                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        dat_Nav = SymIfft(img_Nav,[1 2]);

        % Extracting relative phase maps of 2D navigator images via B0-ESPIRiT
        nk_nav = 7; % empirical kernel size for 2D navigator phase
        img_Nav = zeros(size(img_sense));
        cali_val_set = cell(para.NPhs,para.NSli,para.NRep);
        id_trun_set = cell(para.NPhs,para.NSli,para.NRep);
        for idx_4 = 1:size(img_Nav,4)
            for idx_5 = 1:size(img_Nav,5)
                for idx_6 = 1:size(img_Nav,6)
                    for idx_7 = 1:size(img_Nav,7)
                        for idx_8 = 1:size(img_Nav,8)
                            for idx_9 = 1:size(img_Nav,9)
                                for idx_10 = 1:size(img_Nav,10)
                                    for idx_11 = 1:size(img_Nav,11)
                                        if idx_5 ==1
                                            thres_k_nav = 0.0267;
                                        else
                                            thres_k_nav = round(nk_nav*nk_nav*para.NSeg*0.5);
                                        end

                                        lowRes_dat = dat_Nav(:,:,:,idx_4,idx_5,...
                                            idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);
                                        [eigen_maps,eigen_vals,cali_val,id_trun] = ESPIRiT_map_eigval(lowRes_dat,nk_nav,size(img_sense,1),size(img_sense,2),thres_k_nav);

                                        img_Nav(:,:,:,idx_4,idx_5,...
                                            idx_6,idx_7,idx_8,idx_9,idx_10,idx_11) = eigen_maps(:,:,:,1);

                                        cali_val_set{idx_5,idx_4} = cali_val;% warn: idx_6 and more counters are not considered here
                                        id_trun_set{idx_5,idx_4} = id_trun;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        img_Nav = zeros(size(img_sense));
    end

    if para.bNav ==1 % 2D navigator correction
        if para.coeff_kai ~=0
            apo_MTF = K_2D;
            img_sense = SymFft(SymIfft(img_sense,[1 2]).*apo_MTF,[1 2]);
            apo_MTF = sum(apo_MTF,3);
        end
        if para.bAmpCorr ~=1
            if para.bConjNav ~= 1
                img_sense_temp = img_sense.*exp(1i.*(angle(img_Nav)));
            else
                img_sense_temp = img_sense.*exp(-1i.*(angle(img_Nav)));
            end
        else
            if para.bConjNav ~= 1
                img_sense_temp = img_sense.*img_Nav;
            else
                img_sense_temp = img_sense.*conj(img_Nav);
            end
        end
        NCha=1;
    else
        NCha=1;
        img_sense_temp=img_sense;
    end
    % by the end of the if loop above, keep dim. of img_sense: [Read,Phase,Seg,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]

    % Reorganize the data
    img_sense_temp = permute(img_sense_temp,[1,2,3,11,4,5,6,7,8,9,10]); % [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    dat_sense = SymIfft(img_sense_temp,[1 2]); clear img_sense_temp; % clear to reduce memory

    dat_sense = dat_sense((n_mar+1):(end-n_mar),(n_mar+1):(end-n_mar),:,:,:,:,:,:,:,:,:);

    % eliminate the zero k-space locations, data consistency
    if para.nPFT~=1
        dat_sense(:,1:(loc_pe-1),:,:,:,:,:,:,:,:,:)=0;
    end
    dat_mb = dat_sense;
    img_epi = SymFft(dat_mb,[1 2]);
elseif para.nPFT~=1
    dat_mb = padarray(dat_mb,[0 sum(P==0) 0 0 0 0 0 0 0 0 0],0,'pre');
    img_epi = SymFft(dat_mb,[1 2]);
end

disp('SENSE recon done')

%% Extract inter-segments phase fluctuations

if (para.NSeg>1)  && (para.bIntSCorr==1)

    ind_ol = find(ref_MTF == 2);
    diff_ind_ol = diff(ind_ol);
    breakPoints = [0;find(diff_ind_ol ~= 1);length(ind_ol)];
    loc_ol = cell(1,length(breakPoints)-1);

    for i = 1:(length(breakPoints)-1)
        loc_ol{i} =  ind_ol((breakPoints(i)+1):breakPoints(i+1));
    end

    % crop the overlapped bands
    dat_ol = cell(para.NSli,para.NPhs,para.NSeg,para.NRep);%rtian-4shot
    for idx_rep = 1:para.NRep
        for idx_sl = 1:para.NSli
            for idx_cyc = 1:para.NPhs
                for idx_band = 1:(para.NSeg)%rtian-4shot

                    % dat_mb: [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
                    if para.nPFT ==1
                        switch idx_band % rtian-4shot
                            case 1
                                ind_RO_temp = sum(ref_MTF_shots(:,[1,2]),2)==2;
                                ind_PE_temp = ref_PE_MTF==1;
                                dat_ol{idx_sl,idx_cyc,idx_band,idx_rep} = dat_mb(ind_RO_temp,ind_PE_temp,[1,2],:,idx_sl,idx_cyc,:,:,:,idx_rep,:);
                            case 2
                                ind_RO_temp = sum(ref_MTF_shots(:,[1,3]),2)==2;
                                ind_temp = find(ref_PE_MTF==1);
                                ind_PE_temp = 1:ind_temp(end);
                                dat_ol{idx_sl,idx_cyc,idx_band,idx_rep} = dat_mb(ind_RO_temp,ind_PE_temp,[1,3],:,idx_sl,idx_cyc,:,:,:,idx_rep,:);
                            case 3
                                ind_RO_temp = sum(ref_MTF_shots(:,[3,4]),2)==2;
                                ind_PE_temp = ref_PE_MTF==1;
                                dat_ol{idx_sl,idx_cyc,idx_band,idx_rep} = dat_mb(ind_RO_temp,ind_PE_temp,[3,4],:,idx_sl,idx_cyc,:,:,:,idx_rep,:);
                            case 4
                                ind_RO_temp = sum(ref_MTF_shots(:,[2,4]),2)==2;
                                ind_temp = find(ref_PE_MTF==1);
                                ind_PE_temp = ind_temp(1):length(ref_PE_MTF);
                                dat_ol{idx_sl,idx_cyc,idx_band,idx_rep} = dat_mb(ind_RO_temp,ind_PE_temp,[2,4],:,idx_sl,idx_cyc,:,:,:,idx_rep,:);
                            otherwise
                                error('Error in phase fluctuation estimations: only 4-shot are used')
                        end
                    else
                        dat_ol{idx_sl,idx_cyc,idx_band,idx_rep} = dat_mb(ref_MTF_shots(:,[1 2]),loc_pe:end,[idx_band,idx_band+1],:,idx_sl,idx_cyc,:,:,:,idx_rep,:);
                    end
                end
            end
        end
    end

    % Perform ESPIRiT phase extraction, rtian-4shot
    phase_map = zeros(para.ny_tot,para.nx_tot,2,para.NSli,para.NPhs,para.NSeg,para.NRep);
    cali_val_set = cell((para.NSeg),para.NPhs,para.NSli,para.NRep);
    id_trun_set = cell((para.NSeg),para.NPhs,para.NSli,para.NRep);
    nband_temp = (para.NSeg); % the parfor loop requires iteration within it to be a already calculated positive integer
    for idx_rep = 1:para.NRep
        for idx_sl = 1:para.NSli
            for idx_cyc = 1:para.NPhs
                for idx_band = 1:nband_temp
                    if idx_cyc ==1
                        thres_k = para.thres_k_b0;
                    else
                        thres_k = para.thres_k_diff;
                    end
                    dat_temp = dat_ol{idx_sl,idx_cyc,idx_band,idx_rep};
                    dat_temp = permute(dat_temp,[1,2,4,3]);

                    if para.type_subspace==1
                        [eigen_B0maps,eigen_B0vals,cali_val,id_trun] = ESPIRiT_B0map(dat_temp,para.nky_subspace,para.nkx_subspace,para.ny_tot,para.nx_tot,thres_k);
                        cali_val_set{idx_band,idx_cyc,idx_sl,idx_rep} = cali_val;
                        id_trun_set{idx_band,idx_cyc,idx_sl,idx_rep} = id_trun;
                    elseif para.type_subspace==2
                        [eigen_B0maps,cali_val,id_trun] = GRAPPA_B0map(dat_temp,para.nky_subspace,para.nkx_subspace,para.ny_tot,para.nx_tot,thres_k);
                        cali_val_set{idx_band,idx_cyc,idx_sl,idx_rep} = cali_val;
                        id_trun_set{idx_band,idx_cyc,idx_sl,idx_rep} = id_trun;
                    else
                        error('subspace algorithm types unspecified!')
                    end

                    phase_map(:,:,:,idx_sl,idx_cyc,idx_band,idx_rep) = eigen_B0maps(:,:,:,1);
                end
                disp(['The subspace estimation is done for diffusion direction: ', num2str(idx_cyc), ', slice: ', num2str(idx_sl),', use parfor to speed up.'])
            end
        end
    end

    % calibrate phase maps in different shots
    phase_map(:,:,:,:,:,1,:) = phase_map(:,:,:,:,:,1,:)./phase_map(:,:,1,:,:,1,:).*phase_map(:,:,1,:,:,2,:);
    phase_map(:,:,:,:,:,3,:) = phase_map(:,:,:,:,:,3,:)./phase_map(:,:,1,:,:,3,:).*phase_map(:,:,2,:,:,2,:);
    phase_map(:,:,:,:,:,4,:) = phase_map(:,:,:,:,:,4,:)./phase_map(:,:,1,:,:,4,:).*phase_map(:,:,2,:,:,1,:);

    phase_map_temp = zeros(para.ny_tot,para.nx_tot,para.NSeg,para.NSli,para.NPhs,para.NRep);
    phase_map_temp(:,:,1,:,:,:) = phase_map(:,:,1,:,:,1,:);
    phase_map_temp(:,:,2,:,:,:) = phase_map(:,:,2,:,:,1,:);
    phase_map_temp(:,:,3,:,:,:) = phase_map(:,:,2,:,:,2,:);
    phase_map_temp(:,:,4,:,:,:) = (phase_map(:,:,2,:,:,3,:)+phase_map(:,:,2,:,:,4,:))./2;
    phase_map = phase_map_temp;

    % Take phase or phase+amplitude (not well-tested).
    % phase_map = exp(1i.*angle(phase_map));
    if  para.bAmpCorr ~=1
        phase_map = phase_map./abs(phase_map);
    end
    disp('kernel extraction for B0 fluctuation done.')

    % short axis apodization
    if para.coeff_kai ~=0
        apo_MTF = K_2D;
        img_epi = SymFft(dat_mb.*apo_MTF,[1 2]);
        apo_MTF = sum(apo_MTF,3);
    else
        img_epi = SymFft(dat_mb.*M,[1 2]);
    end

    if para.bMolComb ==0
        % apply phase correction
        img_epi = img_epi.*conj(permute(phase_map,[1,2,3,7,4,5,8,9,10,6,11]));
    else
        phase_map = permute(phase_map,[1,2,3,7,4,5,8,9,10,6,11]);
    end
else
    dat_mb = dat_mb.*M;
    img_epi = SymFft(dat_mb,[1 2]);
end

%% combine shot images

if para.bMolComb ==0

    % sum different shots
    img_epi_comb = sum(img_epi,3);

    % transform to k-space
    dat_epi_comb = SymIfft(img_epi_comb,[1 2]);

    if para.NSeg>1
        if para.bIntSCorr ==1
            % eliminate k-space holes and apply inverse filter
            if para.coeff_kai ~=0
                apo_MTF(apo_MTF<1e-1)=1;
                ref_MTF(ref_MTF<1e-4)=1;
                dat_epi_comb= dat_epi_comb./apo_MTF;
            else
                M_sum(M_sum<1e-4)=1;
                dat_epi_comb= dat_epi_comb./M_sum;
            end
        else
            if para.coeff_kai ~=0
                apo_MTF(apo_MTF<1e-1)=1;
                ref_MTF(ref_MTF<1e-4)=1;
                dat_epi_comb= dat_epi_comb./apo_MTF;
            else
                M_sum(M_sum<1e-4)=1;
                dat_epi_comb= dat_epi_comb./M_sum;
            end
        end
    end

else % solve forward model: dat_epi_comb = E(x), codes can be simplified

    % transform to k-space
    dat_mb = dat_epi_corr;

    % pre-whitening
    if para.bPwhite==1
        dat_mb = permute(dat_mb,[4,1,2,3,5,6,7,8,9,10,11]);
        dat_mb = pagemtimes(L_inv,dat_mb);
        dat_mb = permute(dat_mb,[2,3,4,1,5,6,7,8,9,10,11]);

        load(para.mapS);
        maps_sense = permute(maps_sense,[3,1,2,4]);
        maps_sense = pagemtimes(L_inv,maps_sense);
        maps_sense = permute(maps_sense,[2,3,1,4]);

    end

    dat_mb = padarray(dat_mb,[n_mar 0 0 0 0 0 0 0 0 0 0],0,'both');
    dat_mb = dat_mb.*reshape(ref_MTF_shots,[],1,size(ref_MTF_shots,2));

    ref_MTF_PE_temp = zeros(1,length(P1));
    ref_MTF_PE_temp(1:ind_MTF(end)) = 1;
    M1 = ref_MTF_shots(:,1)*ref_MTF_PE_temp;

    ref_MTF_PE_temp = zeros(1,length(P1));
    ref_MTF_PE_temp(ind_MTF(1):end) = 1;
    M2 = ref_MTF_shots(:,2)*ref_MTF_PE_temp;

    ref_MTF_PE_temp = zeros(1,length(P1));
    ref_MTF_PE_temp(1:ind_MTF(end)) = 1;
    M3 = ref_MTF_shots(:,3)*ref_MTF_PE_temp;

    ref_MTF_PE_temp = zeros(1,length(P1));
    ref_MTF_PE_temp(ind_MTF(1):end) = 1;
    M4 = ref_MTF_shots(:,4)*ref_MTF_PE_temp;

    M=cat(3,M1,M2,M3,M4);

    % Construct Fourier matrix
    E = dftmtx(para.nx_tot+n_mar*2);
    E = circshift(conj(E),[floor(para.nx_tot/2+n_mar) floor(para.nx_tot/2+n_mar)]);
    E = E((n_mar+1):(end-n_mar),:);

    % Eliminate undersampled data
    E_P1 = E(P1==1,:);
    E_P2 = E(P2==1,:);E_P2 = flip(E_P2,1);

    % empirical
    E_P1_temp = E_P1;
    E_P1 = E_P2;
    E_P2=E_P1_temp;
    clear E_P1_temp;

    % dat_mb: [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    % img_sense: [Read,Phase,Seg,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    img_sense = zeros(ny,nx,1,size(dat_mb,5),...
        size(dat_mb,6),size(dat_mb,7),size(dat_mb,8),...
        size(dat_mb,9),size(dat_mb,10),size(dat_mb,11));

    n3_temp  = size(dat_mb,3);
    n5_temp  = size(dat_mb,5);
    n6_temp  = size(dat_mb,6);
    n7_temp  = size(dat_mb,7);
    n8_temp  = size(dat_mb,8);
    n9_temp  = size(dat_mb,9);
    n10_temp  = size(dat_mb,10);
    n11_temp  = size(dat_mb,11);

    if (para.NSeg>1)  && (para.bIntSCorr==1)
        phase_map_temp = phase_map;
    elseif (para.bNav == 1)
        phase_map_temp = img_Nav./(abs(img_Nav)+1e-10);
        if para.bConjNav ==1
            phase_map_temp = conj(phase_map_temp);
        end
        phase_map_temp = permute(phase_map_temp,[1,2,3,7,4,5,8,9,10,6,11]);
    else
        phase_map_temp=1;
    end

    if para.bDisCorr==1
        B0_fmaps_temp = B0_fmaps;
    else
        B0_fmaps_temp=1;
    end

    for idx_5 = 1:n5_temp
        for idx_6 = 1:n6_temp
            for idx_7 = 1:n7_temp
                for idx_8 = 1:n8_temp
                    for idx_9 = 1:n9_temp
                        for idx_10 = 1:n10_temp
                            for idx_11 = 1:n11_temp

                                dat_temp = dat_mb(:,:,:,:,idx_5,...
                                    idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);
                                sz_dat_temp = size(dat_temp);
                                dat_temp = dat_temp(:);

                                temp_sense = reshape(maps_sense(:,:,:,idx_5),ny,nx,1,[]);
                                sz_temp_sense = size(temp_sense);

                                if (para.NSeg>1)  && (para.bIntSCorr==1)
                                    phase_temp = phase_map_temp(:,:,:,:,idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);
                                elseif (para.bNav == 1)
                                    phase_temp = phase_map_temp(:,:,:,:,idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);
                                else
                                    phase_temp=phase_map_temp;
                                end

                                sz_phase_temp = size(phase_temp);

                                if para.bDisCorr==1
                                    B0maps_temp = reshape(B0_fmaps_temp(:,:,idx_5),1,size(B0_fmaps_temp,1),size(B0_fmaps_temp,2));
                                    B0_scale = (-para.ind_kc+(1:sum(P1==1)))./para.deltaTE.*para.ESP;B0_scale = reshape(B0_scale,[],1);
                                    B0maps_temp = exp(1i.*angle(B0maps_temp).*B0_scale);
                                else
                                    B0maps_temp=B0_fmaps_temp;
                                end

                                Afun_reco_4sEPI = @(x,flag) reco_4sEPI(x,flag,temp_sense,phase_temp,B0maps_temp,...
                                    ny,nx,sz_dat_temp,sz_temp_sense,sz_phase_temp,M,ref_MTF_shots,E_P1,E_P2);

                                [x_temp,flag,relres,iter,resvec,lsvec] = lsqr(Afun_reco_4sEPI,dat_temp,[],para.iterLSQR_combShot);
                                % x_temp = pinv(E_sense)*dat_temp;

                                img_sense(:,:,:,idx_5,...
                                    idx_6,idx_7,idx_8,idx_9,idx_10,idx_11) = reshape(x_temp,ny,nx);
                            end
                        end
                    end
                end
            end

        end
    end
    dat_epi_comb = SymIfft(img_sense,[1 2]);
end

%% Apply partial Fourier transform

if para.nPFT ~= 1

    nIter = 10;

    dat_sym = dat_epi_comb;
    dat_sym(:,1:(loc_pe-1),:,:,:,:) = 0;
    dat_sym(:,(end-loc_pe+2):end,:,:,:,:) = 0;

    img_sym = SymFft(dat_sym,[1 2]);
    rho = exp(1i.*angle(img_sym));

    dat_reco = zeros(size(dat_epi_comb));
    for idx_it = 1:nIter

        if idx_it ==1
            dat_reco(:,1:(loc_pe-1),:,:,:,:) = 0;
        else
            % replacing the non sampled region
            dat_reco(:,1:(loc_pe-1),:,:,:,:) = dat_const(:,1:(loc_pe-1),:,:,:,:);
        end
        % enforcing the sampled region
        dat_reco(:,loc_pe:end,:,:,:,:) = dat_epi_comb(:,loc_pe:end,:,:,:,:);

        img_reco = SymFft(dat_reco,[1 2]);
        img_reco = abs(img_reco);
        img_reco = img_reco.*rho;

        dat_const = SymIfft(img_reco,[1 2]);

    end

    dat_reco(:,1:(loc_pe-1),:,:,:,:) = dat_const(:,1:(loc_pe-1),:,:,:,:);
    dat_reco(:,loc_pe:end,:,:,:,:) = dat_epi_comb(:,loc_pe:end,:,:,:,:);

else
    dat_reco = dat_epi_comb;
end

disp('Partial Fourier algorithm done')

%% To obtain the sos images

img_epi_comb = SymFft(dat_reco,[1 2]);

% combine RF channels
if para.bMolComb ==0
    img_epi_ms_sos = squeeze(sqrt(sum(img_epi_comb.*conj(img_epi_comb),4)));
else
    img_epi_ms_sos=squeeze(abs(img_epi_comb));
end

figure
subplot(1,2,1)
imagesc(rot90(sqrt(sum(img_epi_ms_sos(:,:,2,1).*conj(img_epi_ms_sos(:,:,2,1)),4)),-1))
colormap(gray)
title('b = 0 s/mm2')

subplot(1,2,2)
imagesc(rot90(sqrt(sum(img_epi_ms_sos(:,:,2,2:end).*conj(img_epi_ms_sos(:,:,2,2:end)),4)),-1))
colormap(gray)
title('b = 1000 s/mm2, sum, 2 DW directions')

%% End of the script
 
toc(tst_tot);

%% Function to combine 4-shot mosaic segments

function y = reco_4sEPI(x,flag,temp_sense,phase_temp,B0maps_temp,...
    ny,nx,sz_dat_temp,sz_temp_sense,sz_phase_temp,M,ref_MTF_shots,E_P1,E_P2)
if strcmp(flag,'notransp') % Compute A*x
    nband =4;
    % x->2D (ny.nx)
    x = reshape(x,ny,nx);

    % x->image-space phase and sense modulation (ny,nx,nband,nc)
    x = x.*phase_temp.*temp_sense;
    if size(x,3)~=nband
        x=repmat(x,1,1,nband,1);
    end

    % Prepare encoding matrix, DFT phase and B0 off-resonance
    % x-> (nx,nband,nc,y)
    x = permute(x,[2 3 4 1]);

    % y initiatation (kx,nband,nc,ky)
    y = zeros(size(E_P1,1),size(x,2),size(x,3),size(x,4));
    for idx_band = 1:size(x,2)

        % which phase encode direction
        if sum(ismember([1,3],idx_band))
            E_PE_temp = E_P1;
        else
            E_PE_temp = E_P2;
        end

        if sum(size(B0maps_temp)~=1) % distortion correction
            parfor idx_RO = 1:ny
                E_PE = squeeze(B0maps_temp(:,idx_RO,:)).*E_PE_temp;
                y(:,idx_band,:,idx_RO) = E_PE*squeeze(x(:,idx_band,:,idx_RO));
            end
        else
            x_temp = squeeze(x(:,idx_band,:,:));
            x_temp = reshape(x_temp,size(x_temp,1),[]);
            y_temp = E_PE_temp*x_temp;
            y_temp = reshape(y_temp,size(y_temp,1),1,size(temp_sense,4),ny);
            y(:,idx_band,:,:) = y_temp;
        end
    end

    % ->y(y,kx,nband,nc)
    y = permute(y,[4,1,2,3]);

    % y-> (ky,nx,nband,nc)
    y = SymIfft(y,1);

    % data consistency along readout
    y = y.*reshape(ref_MTF_shots,[],1,size(ref_MTF_shots,2));

    % reshape into 1D vector
    y = y(:);

elseif strcmp(flag,'transp') % Compute A'*x

    nPE = size(E_P1,1);

    % hard code
    nband = 4;

    % x(ky,nPE,nband,nc)
    x = reshape(x,ny,nPE,nband,size(temp_sense,4));

    % data consistency along readout
    x = x.*reshape(ref_MTF_shots,[],1,size(ref_MTF_shots,2));

    % ->x(y,nPE,nband,nC)
    x = SymFft(x,1);

    % permute -> y(nPE,nband,nc,ky)
    x = permute(x,[2,3,4,1]);

    % x initialization (x,nband,nc,ky)
    y = zeros(nx,size(x,2),size(x,3),ny);
    for idx_band = 1:size(x,2)

        % which phase encode direction
        if sum(ismember([1,3],idx_band))
            E_PE_temp = E_P1;
        else
            E_PE_temp = E_P2;
        end

        if sum(size(B0maps_temp)~=1) % distortion correction
            parfor idx_RO = 1:ny
                E_PE = squeeze(B0maps_temp(:,idx_RO,:)).*E_PE_temp;
                y(:,idx_band,:,idx_RO) = (E_PE')*squeeze(x(:,idx_band,:,idx_RO));
            end

        else
            x_temp = squeeze(x(:,idx_band,:,:));
            x_temp = reshape(x_temp,size(x_temp,1),[]);
            y_temp = E_PE_temp'*x_temp;
            y_temp = reshape(y_temp,size(y_temp,1),1,size(temp_sense,4),ny);
            y(:,idx_band,:,:) = y_temp;
        end
    end

    % y->(y,x,nband,nc)
    y = permute(y,[4,1,2,3]);

    y = y.*conj(phase_temp).*conj(temp_sense);
    y = mean(y,[3,4]);

    % remove MTF
    y = SymIfft(y,[1 2]);
    M_temp = sum(M,3);M_temp(M_temp<1e-4)=1e6;
    y = y./M_temp;
    y = SymFft(y,[1 2]);

    y=y(:);

end
end


%% Function EPI gridding

function [dat_epi_rg,P,mapS] = gridEPI(dat_epi,k_acq,nOS,nPI,bType,para)

ny_seg = size(k_acq,2);
nRF = size(dat_epi,2);

% The entire k-space range for phase encoding dim.
kxmin = min(k_acq(2,:,:,:,:,:,:,:,:,:,:),[],'all');
kxmax = max(k_acq(2,:,:,:,:,:,:,:,:,:,:),[],'all');

if nPI~=0
    % Estimate phase encoding full grid number
    kx_Step = abs(k_acq(2,1,round(size(k_acq,3)/2),1,1,1,1,1,1,1,1)-k_acq(2,1,round(size(k_acq,3)/2)-1,1,1,1,1,1,1,1,1))/nPI;
    kxmm = max(abs(kxmin),kxmax);
    para.nx_tot = length((-kxmm:kx_Step:kxmm));
    if mod(para.nx_tot,2)==1 % experimental trial, assuming the phaes enocoding steps are designed as even
        para.nx_tot=para.nx_tot-1;
    end

    if para.nx_tot<size(k_acq,3)
        para.nx_tot = size(k_acq,3);
    end

    % MZ manner to calculate cont. k-space window
    kxmax1=kxmax/(para.nx_tot/2-1)*(para.nx_tot/2);
    kx_maxabs=max(kxmax1, -kxmin);
    if bType==1 % acquisition
        kxx= ((-para.nx_tot/2):(para.nx_tot/2-1))/(para.nx_tot/2)*kx_maxabs; % full kx-sample positions
    else % 2D navigator
        kxx= ((-para.nx_tot/2):(para.nx_tot/2))/(para.nx_tot/2)*kx_maxabs; % full kx-sample positions
    end

    % Calculate the sampling pattern
    kx_exp = kxx(1:nPI:end);kx_exp=kx_exp(1:size(k_acq,3)); % actual k-space sampling locations
    kxx = flip(kxx,2);
    kx_exp = flip(kx_exp,2);
    P = ismember(kxx, kx_exp); % The sampling index
else
    P = [];
    para.nx_tot = length(P);
end
% ===============================================================
% Read grid
% ===============================================================
% The entire k-space range for all readout segments
kymin=min(k_acq(1,:,:,:,:,:,:,:,:,:,:),[],'all');
kymax=max(k_acq(1,:,:,:,:,:,:,:,:,:,:),[],'all');

% obtain the k-space step size and estimate the total number of grid
ky_Step = k_acq(1,round(ny_seg/2),1,1,1,1,1,1,1,1,1)-k_acq(1,round(ny_seg/2)-nOS,1,1,1,1,1,1,1,1,1);
para.ny_tot = length((kymin:abs(ky_Step):kymax));
if bType == 1 % acquisition
    mapS = ['maps_Sense_',num2str(para.ny_tot),'_',num2str(para.nx_tot),'_MID',num2str(para.nEPI),'.mat'];
elseif bType ==2 % navigator
    mapS = ['maps_NavSense_',num2str(para.ny_tot),'_',num2str(para.nx_tot),'_MID',num2str(para.nEPI),'.mat'];
end

% MZ manner to calculate con. k-space window
kymax1=kymax/(para.ny_tot/2-1)*(para.ny_tot/2); % this compensates for the non-symmetric center definition in FFT
ky_maxabs=max(kymax1, -kymin);
kyy= ((-para.ny_tot/2):(para.ny_tot/2-1))/(para.ny_tot/2)*ky_maxabs; % kx-sample positions

% Gridding
dat_epi_rg = zeros(para.ny_tot,nRF,...
    size(dat_epi,3),size(dat_epi,4),size(dat_epi,5),...
    size(dat_epi,6),size(dat_epi,7),size(dat_epi,8),...
    size(dat_epi,9),size(dat_epi,10),size(dat_epi,11));

for idx_11 = 1:size(dat_epi,11)
    for idx_10 = 1:size(dat_epi,10)
        for idx_9 = 1:size(dat_epi,9)
            for idx_8 = 1:size(dat_epi,8)
                for idx_7 = 1:size(dat_epi,7)
                    for idx_6 = 1:size(dat_epi,6)
                        for idx_5 = 1:size(dat_epi,5)
                            for idx_4 = 1:size(dat_epi,4)
                                for idx_3 = 1:size(dat_epi,3)
                                    parfor idx_2 = 1:nRF
                                        % Grid, and correct reverse direction, set outside zero
                                        dat_epi_rg(:,idx_2,idx_3,idx_4,idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11) = ...
                                            interp1(k_acq(1,:,idx_3,idx_4,idx_5,idx_6,idx_7,idx_8,1,idx_10,idx_11),... % index 9 is repetition in protocol card, always 1
                                            dat_epi(:,idx_2,idx_3,idx_4,idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11),...
                                            kyy,'spline',0);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

end

%% Function EPI gridding, for interpolatinig on a 2D grid, adapted for phase segments

function [dat_epi_rg,P1,P2,mapS] = grid2DEPI(dat_epi,k_acq,nOS,nPI,bType,para)

ny_seg = size(k_acq,2);
nRF = size(dat_epi,2);

% The entire k-space range for phase encoding dim.
kxmin = min(k_acq(2,:,:,:,:,:,:,:,:,:,:),[],'all');
kxmax = max(k_acq(2,:,:,:,:,:,:,:,:,:,:),[],'all');

if nPI~=0
    % Estimate phase encoding full grid number
    kx_Step = abs(k_acq(2,1,round(size(k_acq,3)/2),1,1,1,1,1,1,1,1)-k_acq(2,1,round(size(k_acq,3)/2)-1,1,1,1,1,1,1,1,1))/nPI;
    % if mod(para.nx_tot,2)==1 % experimental trial, assuming the phaes enocoding steps are designed as even
    %     para.nx_tot=para.nx_tot-1; % commented by rtian-2shot
    % end

    % Define tolerance
    tol_k = kx_Step.*1e-3;

    % Normalize the vectors by rounding to the nearest tolerance
    normalized_vec1 = reshape(round(squeeze(k_acq(2,1,:,1,1,1,1,1,1,1,1)) / tol_k) * tol_k,1,[]);
    normalized_vec2 = reshape(round(squeeze(k_acq(2,1,:,1,1,1,1,1,1,1,2)) / tol_k) * tol_k,1,[]);

    kxx = union(normalized_vec1,normalized_vec2);

    kxx_temp = 1e10.*ones(1, length(kxx) + (length(kxx) - 1) * (nPI-1));
    kxx_temp(1:(nPI):end) = kxx;
    kxx=kxx_temp;clear kxx_temp;

    P1 = ismember(kxx,normalized_vec1);
    P2 = ismember(kxx,normalized_vec2);
    para.nx_tot = length(kxx);
    if para.nx_tot<size(k_acq,3)
        para.nx_tot = size(k_acq,3);
    end
else
    P = [];
    para.nx_tot = length(P);
end
% ===============================================================
% Read grid
% ===============================================================
% The entire k-space range for all readout segments
kymin=min(k_acq(1,:,:,:,:,:,:,:,:,:,:),[],'all');
kymax=max(k_acq(1,:,:,:,:,:,:,:,:,:,:),[],'all');

% obtain the k-space step size and estimate the total number of grid
ky_Step = k_acq(1,round(ny_seg/2),1,1,1,1,1,1,1,1,1)-k_acq(1,round(ny_seg/2)-nOS,1,1,1,1,1,1,1,1,1);
para.ny_tot = length((kymin:abs(ky_Step):kymax));
if bType == 1 % acquisition
    mapS = ['maps_Sense_',num2str(para.ny_tot),'_',num2str(para.nx_tot),'_MID',num2str(para.nEPI),'.mat'];
elseif bType ==2 % navigator
    mapS = ['maps_NavSense_',num2str(para.ny_tot),'_',num2str(para.nx_tot),'_MID',num2str(para.nEPI),'.mat'];
end

% MZ manner to calculate con. k-space window
kymax1=kymax/(para.ny_tot/2-1)*(para.ny_tot/2); % this compensates for the non-symmetric center definition in FFT
ky_maxabs=max(kymax1, -kymin);
kyy= ((-para.ny_tot/2):(para.ny_tot/2-1))/(para.ny_tot/2)*ky_maxabs; % kx-sample positions

% Gridding
if sum(P1)~=sum(P2)
    error('error:phase encoded step numer in 2 bands are not identical!')
end
dat_epi_rg = zeros(para.ny_tot,nRF,...
    sum(P1),size(dat_epi,4),size(dat_epi,5),...
    size(dat_epi,6),size(dat_epi,7),size(dat_epi,8),...
    size(dat_epi,9),size(dat_epi,10),size(dat_epi,11));

for idx_11 = 1:size(dat_epi,11)
    for idx_10 = 1:size(dat_epi,10)
        for idx_9 = 1:size(dat_epi,9)
            for idx_8 = 1:size(dat_epi,8)
                for idx_7 = 1:size(dat_epi,7)
                    for idx_6 = 1:size(dat_epi,6)
                        for idx_5 = 1:size(dat_epi,5)
                            for idx_4 = 1:size(dat_epi,4)
                                for idx_3 = 1:size(dat_epi,3)
                                    parfor idx_2 = 1:nRF
                                        % Grid, and correct reverse direction, set outside zero
                                        dat_epi_rg(:,idx_2,idx_3,idx_4,idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11) = ...
                                            interp1(k_acq(1,:,idx_3,idx_4,idx_5,idx_6,idx_7,idx_8,1,idx_10,idx_11),... % index 9 is repetition in protocol card, always 1
                                            dat_epi(:,idx_2,idx_3,idx_4,idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11),...
                                            kyy,'spline',0);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

end

%% Function to apply even/odd echoes correction

function [dat_epi_1Dcorr,img_epi_1Dcorr] = corr_EvenOddEchoes(dat_epi,offset_avg,ref_MTF_shots,ny_tot,ny_echoes,bIntpHR,bSP)

if bIntpHR==1
    offset_avg = SymIfft(offset_avg,1);
    if mod(ny_echoes,2)==1
        offset_avg = padarray(offset_avg,ceil((ny_tot-ny_echoes)/2),0,'pre');
        offset_avg = padarray(offset_avg,floor((ny_tot-ny_echoes)/2),0,'post');
    else % mod(ny_echoes,2)==0
        offset_avg = padarray(offset_avg,ceil((ny_tot-ny_echoes)/2),0,'post');
        offset_avg = padarray(offset_avg,floor((ny_tot-ny_echoes)/2),0,'pre');
    end
    offset_avg = SymFft(offset_avg,1);
    dat_epi_1Dcorr = SymFft(dat_epi,1); %[read,phase,seg,rf,slice,cycle]
    if bSP==0
        dat_epi_1Dcorr(:,2:2:end,:,:,:,:,:,:,:,:,:) = dat_epi_1Dcorr(:,2:2:end,:,:,:,:,:,:,:,:,:).*permute(exp(1i.*angle(offset_avg)),[1 3 5 2 4]);
    elseif bSP==3 %rtian-4shot, correct even echoes for band 1-2, odd echoes for band 3-4
        dat_epi_1Dcorr(:,2:2:end,[1,2],:,:,:,:,:,:,:,:) = dat_epi_1Dcorr(:,2:2:end,[1,2],:,:,:,:,:,:,:,:).*permute(exp(1i.*angle(offset_avg)),[1 3 5 2 4]);
        dat_epi_1Dcorr(:,1:2:end,[3,4],:,:,:,:,:,:,:,:) = dat_epi_1Dcorr(:,1:2:end,[3,4],:,:,:,:,:,:,:,:).*permute(exp(1i.*angle(offset_avg)),[1 3 5 2 4]);
    else
        if mod((size(dat_epi,3)-1)/2,2)==1
            dat_epi_1Dcorr(:,2:2:end,2:2:end,:,:,:,:,:,:,:,:) = dat_epi_1Dcorr(:,2:2:end,2:2:end,:,:,:,:,:,:,:,:).*permute(exp(1i.*angle(offset_avg)),[1 3 5 2 4]);
            dat_epi_1Dcorr(:,1:2:end,1:2:end,:,:,:,:,:,:,:,:) = dat_epi_1Dcorr(:,1:2:end,1:2:end,:,:,:,:,:,:,:,:).*permute(exp(1i.*angle(offset_avg)),[1 3 5 2 4]);
        else
            dat_epi_1Dcorr(:,2:2:end,1:2:end,:,:,:,:,:,:,:,:) = dat_epi_1Dcorr(:,2:2:end,1:2:end,:,:,:,:,:,:,:,:).*permute(exp(1i.*angle(offset_avg)),[1 3 5 2 4]);
            dat_epi_1Dcorr(:,1:2:end,2:2:end,:,:,:,:,:,:,:,:) = dat_epi_1Dcorr(:,1:2:end,2:2:end,:,:,:,:,:,:,:,:).*permute(exp(1i.*angle(offset_avg)),[1 3 5 2 4]);
        end
    end

    img_epi_1Dcorr = SymFft(dat_epi_1Dcorr,2);%[read,phase,seg,rf,slice,cycle]
    dat_epi_1Dcorr = SymIfft(dat_epi_1Dcorr,1);
else % does not work for switched polarity yet
    dat_epi_1Dcorr = dat_epi;
    for idx_pe = 2:2:size(dat_epi,2)
        for idx_shot = 1:size(dat_epi,3)
            for idx_rf = 1:size(dat_epi,4)
                for idx_5 = 1:size(dat_epi,5)
                    for idx_6 = 1:size(dat_epi,6)
                        for idx_7 = 1:size(dat_epi,7)
                            for idx_8 = 1:size(dat_epi,8)
                                for idx_9 = 1:size(dat_epi,9)
                                    for idx_10 = 1:size(dat_epi,10)
                                        for idx_11 = 1:size(dat_epi,11)

                                            dat_epi_temp = dat_epi(:,idx_pe,idx_shot,idx_rf,...
                                                idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);

                                            ref_MTF_temp = ref_MTF_shots(:,idx_shot);
                                            dat_epi_temp = dat_epi_temp(ref_MTF_temp~=0);
                                            dat_epi_temp = SymFft(dat_epi_temp,1); %[read,phase,seg,rf,slice,cycle]

                                            offset_temp = offset_avg;
                                            if size(dat_epi_temp,1)~= size(offset_avg,1)
                                                offset_temp = imresize(offset_temp,[size(dat_epi_temp,1),1]);
                                            end
                                            dat_epi_temp = dat_epi_temp.*permute(exp(1i.*angle(offset_temp)),[1 3 5 2 4]);

                                            temp = dat_epi(:,idx_pe,idx_shot,idx_rf);
                                            temp(ref_MTF_temp~=0) = SymIfft(dat_epi_temp,1);

                                            dat_epi_1Dcorr(:,idx_pe,idx_shot,idx_rf,...
                                                idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11)=temp;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    img_epi_1Dcorr = SymFft(dat_epi_1Dcorr,[1 2]);%[read,phase,seg,rf,slice,cycle]
end

end


%% Function: ESPIRiT-type operator to extract shot-to-shot phase fluctuation maps

function [eigen_maps,eigen_vals,cali_val,id_trun] = ESPIRiT_B0map(lowRes_dat,nky_subspace,nkx_subspace,ny,nx,thres_k)

% Parameters for calibration region
[ny_low, nx_low, nc, nb] = size(lowRes_dat);
ns_y = ny_low-nky_subspace+1;
ns_x = nx_low-nkx_subspace+1;

ny_zp = ceil((ny-nky_subspace)/2);
nx_zp = ceil((nx-nkx_subspace)/2);

y_zp_odd = ceil((ny-nky_subspace)/2) ~= (ny-nky_subspace)/2;
x_zp_odd = ceil((nx-nkx_subspace)/2) ~= (nx-nkx_subspace)/2;

cali_matx = zeros(ns_y*ns_x,nc,nky_subspace*nkx_subspace*nb);

if size(cali_matx,1)<size(cali_matx,3)
    disp('Warning: No enough shifts to calculate subspace map with this kernel size, return.')
    id_trun=0;
    eigen_maps = zeros(ny,nx,nb,nb);
    eigen_vals = zeros(ny,nx,nb);
    cali_val=0;
    return;
end

% Formalize the calibration matrix
for index_band = 1:nb
    x_cali = (1:nky_subspace*nkx_subspace) + (index_band-1)*nky_subspace*nkx_subspace;
    for idx_rf = 1:nc
        for index_shift_x = 1:ns_x
            x_lowRes = (1:nkx_subspace) + index_shift_x-1;
            for index_shift_y = 1:ns_y
                y_cali = index_shift_x + (index_shift_y-1)*ns_x;
                y_lowRes = (1:nky_subspace) + index_shift_y-1;
                dat_temp = lowRes_dat(y_lowRes,x_lowRes,idx_rf,index_band);
                cali_matx(y_cali,idx_rf,x_cali) = reshape(dat_temp,[1 nky_subspace*nkx_subspace]);
            end
        end
    end
end
cali_matx = reshape(cali_matx,[ns_y*ns_x*nc,nky_subspace*nkx_subspace*nb]);

% SVD of the calibration matrix for thresholding
[U_cali,S_cali,V_cali] = svd(cali_matx,'econ');
cali_val = diag(S_cali);
% different range of thres_k for different thresholding
if (thres_k>0) && (thres_k<1)
    id_trun = find(cali_val>=cali_val(1)*thres_k, 1, 'last' );
elseif thres_k>=1
    id_trun = thres_k;
else
    error('subspace truncation index cannot be equal or smaller than 0')
end

% Empirical, always cut at least one singular vector for stability
if id_trun==nky_subspace*nkx_subspace*nb
    id_trun = id_trun-1;
end

% Reshape to k-space filters
cali_filters = reshape(V_cali,[nky_subspace nkx_subspace nb nky_subspace*nkx_subspace*nb]);
cali_filters = cali_filters(:,:,:,1:id_trun);
cali_filters = permute(cali_filters,[1,2,4,3]);
cali_filters = reshape(cali_filters,nky_subspace*nkx_subspace*id_trun,nb);
[U_filter,S_filter,V_filter] = svd(cali_filters,'econ');
cali_filters = cali_filters*V_filter;
cali_filters = reshape(cali_filters,[nky_subspace nkx_subspace id_trun nb]);
cali_filters = permute(cali_filters,[1,2,4,3]);

% convert to image domain
img_kernel = zeros(ny,nx,nb,id_trun);
for index_eigen = 1:id_trun
    cali_filters_temp = cali_filters(:,:,:,index_eigen);
    cali_filters_temp = conj(flip(flip(cali_filters_temp,1),2));
    cali_filters_temp_zp = padarray(cali_filters_temp,[ny_zp nx_zp],0,'both');
    if y_zp_odd == 1
        cali_filters_temp_zp = cali_filters_temp_zp(1:end-1,:,:);
    end
    if x_zp_odd == 1
        cali_filters_temp_zp = cali_filters_temp_zp(:,1:end-1,:);
    end
    img_kernel(:,:,:,index_eigen) = SymIfft(cali_filters_temp_zp,[1 2]);
end

% SVD in the image domain
eigen_maps = zeros(ny,nx,nb,min(id_trun,nb));
eigen_vals = zeros(ny,nx,min(id_trun,nb));
for id_y = 1:ny
    for id_x =  1:nx
        img_kernel_yx = squeeze(img_kernel(id_y,id_x,:,:));
        [U_img,S_img,V_img] = svd(img_kernel_yx,'econ');
        ph_ref_rad = angle(U_img(1,:));
        ph_comp = repmat(exp(-1i*ph_ref_rad),[nb 1]);
        U_img = V_filter*(U_img.*ph_comp);
        S_img = diag(S_img);
        eigen_maps(id_y,id_x,:,:) = U_img;
        eigen_vals(id_y,id_x,:) = S_img;
    end
end
end

%% Function: GRAPPA-type operator to extract shot-to-shot phase fluctuation maps, codes could be simplified

function [phase_maps,cali_val,id_trun] = GRAPPA_B0map(lowRes_dat,nky_subspace,nkx_subspace,ny,nx,thres_k)

% Parameters for calibration region
[ny_low, nx_low, nc, nb] = size(lowRes_dat);
ns_y = ny_low-nky_subspace+1;
ns_x = nx_low-nkx_subspace+1;

ny_zp = ceil((ny-nky_subspace)/2);
nx_zp = ceil((nx-nkx_subspace)/2);

y_zp_odd = ceil((ny-nky_subspace)/2) ~= (ny-nky_subspace)/2;
x_zp_odd = ceil((nx-nkx_subspace)/2) ~= (nx-nkx_subspace)/2;

cali_matx = zeros(ns_y*ns_x,nc,nky_subspace*nkx_subspace*nb);
est_matx = zeros(ns_y*ns_x,nc,nb);

if size(cali_matx,1)<(size(cali_matx,3)/nb)
    disp('Warning: No enough shifts to calculate subspace map with this kernel size, return.')
    id_trun=0;
    eigen_maps = zeros(ny,nx,nb,nb);
    eigen_vals = zeros(ny,nx,nb);
    cali_val=0;
    return;
end

% Construct the calibration matrix
for index_band = 1:nb
    x_cali = (1:nky_subspace*nkx_subspace) + (index_band-1)*nky_subspace*nkx_subspace;
    for idx_rf = 1:nc
        for index_shift_x = 1:ns_x
            x_lowRes = (1:nkx_subspace) + index_shift_x-1;
            for index_shift_y = 1:ns_y
                y_cali = index_shift_x + (index_shift_y-1)*ns_x;
                y_lowRes = (1:nky_subspace) + index_shift_y-1;
                dat_temp = lowRes_dat(y_lowRes,x_lowRes,idx_rf,index_band);
                cali_matx(y_cali,idx_rf,x_cali) = reshape(dat_temp,[1 nky_subspace*nkx_subspace]);

                cen_y = floor(nky_subspace/2)+1;
                cen_x = floor(nkx_subspace/2)+1;
                est_matx(y_cali,idx_rf,index_band) = lowRes_dat(y_lowRes(cen_y),x_lowRes(cen_x),idx_rf,index_band);
            end
        end
    end
end

cali_matx = reshape(cali_matx,[ns_y*ns_x*nc,nky_subspace*nkx_subspace*nb]);
est_matx = reshape(est_matx,[ns_y*ns_x*nc,nb]);

% one-way interpolation
cali_matx_1 = cali_matx(:,1:nky_subspace*nkx_subspace);
est_matx_1 = est_matx(:,2);
if thres_k~=0
    [u,s,v] = svd(cali_matx_1);
    cali_val = diag(s);
    if (thres_k>0) && (thres_k<1)
        id_trun = find(cali_val>=cali_val(1)*thres_k, 1, 'last' );
    elseif thres_k>=1
        id_trun = thres_k;
    else
        error('subspace truncation index cannot be equal or smaller than 0')
    end

    s = pinv(s); s(id_trun:end,id_trun:end)=0;
    inv_kss = v*s*u';
else
    % pseudoinverse, or tsvd without truncations
    % inv_kss = pinv(cali_matx_1);
    [u,s,v] = svd(cali_matx_1);
    cali_val = diag(s);
    id_trun = length(cali_val)+1;
    s = pinv(s);
    inv_kss = v*s*u';
end
coeff_1 = inv_kss*est_matx_1;

% other solvers possible, e.g., lsqr
% [coeff_1,flag,relres,iter,resvec,lsvec] = lsqr(cali_matx_1,est_matx_1);

% the other direction of interpolation
cali_matx_2 = cali_matx(:,(nky_subspace*nkx_subspace+1):end);
est_matx_2 = est_matx(:,1);
if thres_k~=0
    [u,s,v] = svd(cali_matx_2);
    cali_val = diag(s);
    if (thres_k>0) && (thres_k<1)
        id_trun = find(cali_val>=cali_val(1)*thres_k, 1, 'last' );
    elseif thres_k>=1
        id_trun = thres_k;
    else
        error('subspace truncation index cannot be equal or smaller than 0')
    end
    s = pinv(s); s(id_trun:end,id_trun:end)=0;
    inv_kss = v*s*u';
else
    %inv_kss = pinv(cali_matx_2);
    [u,s,v] = svd(cali_matx_2);
    cali_val = diag(s);
    id_trun = length(cali_val)+1;
    s = pinv(s);
    inv_kss = v*s*u';
end
coeff_2 = inv_kss*est_matx_2;
% [coeff_2,flag,relres,iter,resvec,lsvec] = lsqr(cali_matx_2,est_matx_2);

% kernel converted to image-space maps
coeff_1 = reshape(coeff_1,[nky_subspace,nkx_subspace]);
coeff_2 = reshape(coeff_2,[nky_subspace,nkx_subspace]);
coeff_1 = conj(flip(flip(coeff_1,1),2));
coeff_2 = conj(flip(flip(coeff_2,1),2));
coeff_1 = padarray(coeff_1,[ny_zp nx_zp],0,'both');
coeff_2 = padarray(coeff_2,[ny_zp nx_zp],0,'both');
if y_zp_odd == 1
    coeff_1 = coeff_1(1:end-1,:,:);
    coeff_2 = coeff_2(1:end-1,:,:);
end
if x_zp_odd == 1
    coeff_1 = coeff_1(:,1:end-1);
    coeff_2 = coeff_2(:,1:end-1);
end
coeff_1 = SymIfft(coeff_1,[1 2]);
coeff_2 = SymIfft(coeff_2,[1 2]);
coeff_2 = 1./coeff_2;
phase_rel = 0.5.*(real(coeff_1)+real(coeff_2))+0.5.*1i.*(imag(coeff_1)+imag(coeff_2));

phase_maps = cat(3,phase_rel,ones(size(phase_rel)));

end

%% Get the last struct in case the raw data has 1 struct or a cell array with 2 structs

function outputStruct = getLastStruct(inputVar)
if isstruct(inputVar)
    % If inputVar is a single struct, return it
    outputStruct = inputVar;
elseif iscell(inputVar) && all(cellfun(@isstruct, inputVar)) && (length(inputVar) == 2)
    % If inputVar is a cell array of structs, return the last struct
    outputStruct = inputVar{end};
else
    error('Input must be a struct or a cell array of 2 structs.');
end
end


