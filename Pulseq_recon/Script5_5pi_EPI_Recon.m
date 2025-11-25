%% Script description

% This demo script can reconstruct data acquired from the 5-shot phase-interleaved sequences,
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
% 5-shot, phase-interleaved, 2-slice, b=0,1000s/mm2, 2D navigator ON

%% User input: enter sequence/recon parameters

% Experiment number
para.nEPI = 226; % EPI scan MID
para.nNS = 218; % Noise scan MID

% Intershot correction for phase interleave EPI
para.nshot_pred = 5; % shot number of phase interleaves
para.bIntSCorr = 1; % 1: Apply inter-segment correction, 0: does not apply
para.type_subspace = 1; % 1: ESPIRiT 2: MUSE (with Total variation for phase extraction)
para.nky_cali = 50; % readout calibration region for self-navigation
para.nkx_cali = 50; % phase calibration region for self-navigation
para.nky_subspace=7; % kernel length along y
para.nkx_subspace=7; % kernel length along x
para.thres_k_diff = round(para.nky_subspace*para.nkx_subspace*para.nshot_pred*0.5); % empirical kernel threshold, b=0 s/mm2
para.thres_k_b0 = round(para.nky_subspace*para.nkx_subspace*para.nshot_pred*0.5); % empirical kernel threshold, b~=0 s/mm2

% Parallel imaging
para.bPwhite = 1; % prewhitening for imaging data
para.iterLSQR_sepShot = 10; % Max. iteration for LSQR, intermediate PI recon.
para.iterLSQR_combShot = 4; % Max. iteration for LSQR, combining shots

% 2D Navigator， check the necessity of phase conjugate, if for GE-navigator
para.bNav = 0; % apply 2D navigator
para.thres_k_diff_nav = round(para.nky_subspace*para.nkx_subspace*para.nshot_pred*0.5); % empirical kernel threshold

% Common sequence parameter
para.bAmpCorr = 0; % 1: shot-to-shot amplitude correction included. 0: otherwise
para.nOS = 1; % readout oversampling factor
para.bSP = 0; % readout polarity switch
para.nPI = 1; % parallel imaging acceleration factor
para.nPFT = 6/8; % Partial Fourier
para.bDisCorr = 1; % 1: Use distortion correction, else: no distortion correction.
para.deltaTE = 3; % For 2 GRE field mapping scans, TE2-TE1 in ms
para.MID_FieldMap = [220,221]; % The MIDs for GRE scans with 2 echo times.
para.nCali = 219; % GRE calibration scan for parallel imaging
para.kB0 = 20; % kernel window size of B0 mapping. Try smaller values (e.g., 10) for reducing recon instability given too strong distortion.

% Scan dependent sequence parameter
switch para.nEPI
    case 226
        para.seq_file_path = 'rtian_05_5pi_SE_epi_b1k.seq'; % multi-shot EPI sequence file
        nProj = 231; % MID, 1D navigator scan for odd-even echo correction
        para.seq_proj_file_path = 'rtian_05_proj.seq'; % sequence file for 1D navigator scan
        para.ESP = 1.1; % ESP [ms], see printout of the Pulseq sequence
        para.ind_kc = 23; % index of the central k-space phase encoded lines, for distortion corr. in SE-EPI
        para.ind_Nav_kc = 20; % index of the central k-space phase encoded lines, for distortion corr. in SE-EPI
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

if para.bNav==1
    dat_2DNav = file_epi_ms.phasecor();
    dat_2DNav = permute(dat_2DNav,[1,2,3,4,10,6,7,8,9,5,11]);% Switch SET and SLC label, used to avoid pulseq s
    dat_2DNav = dat_2DNav(:,:,:,:,perm_order,:,:,:,:,:,:);
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
if para.bNav==1 % separate navigator data from imaging data
    para.NLin_2DNav = size(k_acq,3)-para.NLin;
    k_2DNav = k_acq(:,:,end-para.NLin_2DNav+1:end,:,:,:,:,:,:,:,:);
end

% imaging trajectory
k_acq = k_acq(:,:,1:para.NLin,:,:,:,:,:,:,:,:);

% EPI gridding
[dat_epi_rg,P,para.mapS] = gridEPI(dat_epi,k_acq,para.nOS,para.nPI*para.NSeg,1,para);
para.ny_tot = size(dat_epi_rg,1);
para.nx_tot = length(P);

% 2D navigator gridding
if para.bNav==1
    [dat_2DNav_rg,P_Nav,para.mapS_Nav] = gridEPI(dat_2DNav,k_2DNav,para.nOS,para.nPI,2,para);
    ny_Nav = size(dat_2DNav_rg,1);
    nx_Nav = length(P_Nav);
end

% k-space coverage
ref_MTF_shots = double(squeeze((dat_epi_rg(:,1,1,1,1,1,1,1,1,1,:)~=0)));
ref_MTF = sum(ref_MTF_shots,2);

dat_epi = dat_epi_rg; clear dat_epi_rg; % clear to reduce memory
dat_epi = permute(dat_epi,[1,3,11,2,5,7,4,6,8,9,10]); %-> update dim. [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]

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

% Steps of phase encoded lines to take into account
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
[dat_epi_corr,img_epi_1Dcorr] = corr_EvenOddEchoes(dat_epi,offset_avg,ref_MTF_shots,para.ny_tot,ny_echoes,bIntpHR,para.bSP);
if para.bNav==1
    [dat_nav_1Dcorr,img_nav_1Dcorr] = corr_EvenOddEchoes(dat_2DNav,offset_avg,ref_MTF_shots,ny_echoes,ny_echoes,bIntpHR,0); % read seg EPI: 2 input of ny_echoes to avoid interpolation
end

% The N/2 ghost corrected image
img_epi_1Dcorr_sos = squeeze(sqrt(sum(img_epi_1Dcorr.*conj(img_epi_1Dcorr),4)));

disp('Naive N/2 ghost correction done.')

%% Pre-process data

% Transform the k-space dat to image space
img_epi = SymFft(dat_epi_corr,[1 2]);

% start from multi segment k-space data
dat_mb = dat_epi_corr;
loc_pe = find(P==1);loc_pe=loc_pe(1);

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
        maps_mask = maps_eigval(:,:,1,:)>(0.9998.*max(maps_eigval(:,:,1,:),[],[1,2]));
        save(para.mapS,'maps_sense','maps_eigval','maps_mask')
    else
        load(para.mapS);
    end

    % Load the GRE images for estimating B0 maps
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
            maps_Nav_mask = maps_Nav_eigval(:,:,1,:)>(0.9998.*max(maps_Nav_eigval(:,:,1,:),[],[1,2]));
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
        E_Nav = E_Nav(P_Nav==1,:);
        E_Nav = repmat(E_Nav,[para.NCha,1]);
    end

    % dat_mb: [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    % img_sense: [Read,Phase,Seg,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    img_sense = zeros(ny,nx,para.NSeg,size(dat_mb_y_kx,5),...
        size(dat_mb_y_kx,6),size(dat_mb_y_kx,7),size(dat_mb_y_kx,8),...
        size(dat_mb_y_kx,9),size(dat_mb_y_kx,10),size(dat_mb_y_kx,11));
    if para.bNav==1
        img_Nav = zeros(ny_Nav,nx_Nav,para.NSeg,size(dat_nav_mb_y_kx,5),...
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
                                        temp_sense = repmat(temp_sense,[sum(P==1),1,1]);
                                        temp_sense = permute(temp_sense,[1,3,2]);
                                        temp_sense = reshape(temp_sense,[],nx);

                                        if para.bDisCorr==1
                                            B0maps_temp = B0_fmaps_temp(idx_RO,:,idx_5);
                                            B0_scale = (-para.ind_kc+(1:sum(P==1)))./para.deltaTE.*para.ESP;B0_scale = reshape(B0_scale,[],1);
                                            B0maps_temp = exp(1i.*angle(B0maps_temp).*B0_scale);
                                            B0maps_temp = repmat(B0maps_temp,[para.NCha,1]);
                                        else
                                            B0maps_temp=B0_fmaps_temp;
                                        end

                                        % Select phase encoded steps for all interleave
                                        E_temp = [];
                                        P_temp = circshift(P==1,-(idx_3-1)*para.nPI);
                                        E_temp = E(P_temp,:);

                                        % Copy for sense
                                        E_temp = repmat(E_temp,[para.NCha,1]);

                                        E_sense = E_temp.*temp_sense.*B0maps_temp;

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
        if size(img_Nav,3)~=size(dat_Nav,3)
            img_Nav = repmat(img_Nav,[1,1,size(dat_Nav,3),1,1,1,1,1,1,1,1]);
        end
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
                                            thres_k_nav = para.thres_k_diff_nav;
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

    % Just inherit the variable. 2D navigator will be applied in the
    % last step of forward model to combine shots
    img_sense_temp = img_sense;
    NCha=1;
    % by the end of the if loop above, keep dim. of img_sense: [Read,Phase,Seg,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]

    % Reorganize the data
    img_sense_temp = permute(img_sense_temp,[1,2,3,11,4,5,6,7,8,9,10]); % [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
    dat_sense = SymIfft(img_sense_temp,[1 2]); clear img_sense_temp; % clear to reduce memory

    dat_sense = dat_sense((n_mar+1):(end-n_mar),(n_mar+1):(end-n_mar),:,:,:,:,:,:,:,:,:);

    % Eliminate the zero k-space locations, data consistency
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

    if ismember(para.type_subspace,[1]) % kernel B0 estimation
        ncen_RO = floor(para.ny_tot/2)+1;
        loc_ol_ro = (ncen_RO-floor(para.nky_cali/2)) : (ncen_RO-1+ceil(para.nky_cali/2));
        ncen_PE = floor(para.nx_tot/2)+1;
        loc_ol_pe = (ncen_PE-floor(para.nkx_cali/2)) : (ncen_PE-1+ceil(para.nkx_cali/2));
        if loc_ol_pe(1)<loc_pe
            loc_ol_pe= loc_pe:loc_ol_pe(end);
        end

        % crop the overlapped bands
        dat_ol = cell(para.NSli,para.NPhs,1,para.NRep);
        for idx_rep = 1:para.NRep
            for idx_sl = 1:para.NSli
                for idx_cyc = 1:para.NPhs
                    % dat_mb: [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
                    dat_ol{idx_sl,idx_cyc,1,idx_rep} = dat_mb(loc_ol_ro,loc_ol_pe,:,:,idx_sl,idx_cyc,:,:,:,idx_rep,:);
                end
            end
        end

        % Perform ESPIRiT phase extraction
        phase_map = zeros(para.ny_tot,para.nx_tot,para.NSeg,para.NSli,para.NPhs,1,para.NRep); % modify for phase interleave EPI
        cali_val_set = cell(1,para.NPhs,para.NSli,para.NRep);
        id_trun_set = cell(1,para.NPhs,para.NSli,para.NRep);
        for idx_rep = 1:para.NRep
            for idx_sl = 1:para.NSli
                for idx_cyc = 1:para.NPhs

                    if idx_cyc ==1
                        thres_k = para.thres_k_b0;
                    else
                        thres_k = para.thres_k_diff;
                    end

                    dat_temp = dat_ol{idx_sl,idx_cyc,1,idx_rep};
                    dat_temp = permute(dat_temp,[1,2,4,3]);

                    if para.type_subspace==1
                        [eigen_B0maps,eigen_B0vals,cali_val,id_trun] = ESPIRiT_B0map(dat_temp,para.nky_subspace,para.nkx_subspace,para.ny_tot,para.nx_tot,thres_k);
                        cali_val_set{1,idx_cyc,idx_sl,idx_rep} = cali_val;
                        id_trun_set{1,idx_cyc,idx_sl,idx_rep} = id_trun;
                    else
                        error('subspace algorithm types unspecified!')
                    end

                    phase_map(:,:,:,idx_sl,idx_cyc,1,idx_rep) = eigen_B0maps(:,:,:,1);

                    % disp(['The subspace estimation is done for diffusion direction: ', num2str(idx_cyc), ', slice: ', num2str(idx_sl)])
                end
            end
        end

    else % total variation B0 estimation, could be optimized further.

        ncen_RO = floor(para.ny_tot/2)+1;
        loc_ol_ro = (ncen_RO-floor(para.nky_cali/2)) : (ncen_RO-1+ceil(para.nky_cali/2));
        ncen_PE = floor(para.nx_tot/2)+1;
        loc_ol_pe = (ncen_PE-floor(para.nkx_cali/2)) : (ncen_PE-1+ceil(para.nkx_cali/2));
        if loc_ol_pe(1)<loc_pe
            loc_ol_pe= loc_pe:loc_ol_pe(end);
        end

        dat_ol = dat_mb(loc_ol_ro,loc_ol_pe,:,:,:,:,:,:,:,:,:);

        coeff_kai_sn = 5;
        kai_nav_x = kaiser(size(dat_ol,1),coeff_kai_sn);
        kai_nav_y = kaiser(size(dat_ol,2),coeff_kai_sn);
        filter_nav = kai_nav_x.*kai_nav_y';
        dat_ol = dat_ol.*filter_nav;

        phase_map = zeros(para.ny_tot,para.nx_tot,para.NSeg,para.NSli,para.NPhs,1,para.NRep); % modify for phase interleave EPI
        img_epi_TV = SymFft(dat_ol,[1 2]);
        para.TV_Theta=16;
        for idx_3=1:size(img_epi_TV,3)
            for idx_5=1:size(img_epi_TV,5)
                for idx_6=1:size(img_epi_TV,6)
                    for idx_10=1:size(img_epi_TV,10)
                        phase_map_temp = ROFdenoise(img_epi_TV(:,:,idx_3,1,idx_5,idx_6,1,1,1,idx_10,1),para.TV_Theta);
                        phase_map(:,:,idx_3,idx_5,idx_6,1,idx_10) = imresize(real(phase_map_temp),[para.ny_tot,para.nx_tot])+...
                            1i.*imresize(imag(phase_map_temp),[para.ny_tot,para.nx_tot]);
                    end
                end
            end
        end

    end

    % Take phase or phase+amplitude (not well-tested).
    % phase_map = exp(1i.*angle(phase_map));
    if  para.bAmpCorr ~=1
        phase_map = phase_map./(abs(phase_map)+1e-10);
    end

    disp('Extraction for B0 fluctuation done.')

    phase_map = permute(phase_map,[1,2,3,4,5,8,9,10,7,6,11]);

end

%% Combine shot images

% Solve forward model: dat_epi_comb = E(x), codes could be simplified.

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
dat_mb_y_kx = SymFft(dat_mb,1);

% Construct Fourier matrix
E = dftmtx(para.nx_tot+n_mar*2);
E = circshift(conj(E),[floor(para.nx_tot/2+n_mar) floor(para.nx_tot/2+n_mar)]);
E = E((n_mar+1):(end-n_mar),:);

% Select phase encoded steps for all interleave
E_temp = [];
for idx_seg = 1:para.NSeg
    P_temp = circshift(P==1,-(idx_seg-1)*para.nPI);
    E_temp = cat(1,E_temp,E(P_temp,:));
end
E=E_temp;clear E_temp;

% Copy for sense
E = repmat(E,[para.NCha,1]);

% dat_mb: [Read,Phase,Seg,RF,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
% img_sense: [Read,Phase,Seg,2D Slice,NPhs(diff),NPar,NAve,NEco,NRep,NSet]
img_sense = zeros(ny,nx,1,size(dat_mb_y_kx,5),...
    size(dat_mb_y_kx,6),size(dat_mb_y_kx,7),size(dat_mb_y_kx,8),...
    size(dat_mb_y_kx,9),size(dat_mb_y_kx,10),size(dat_mb_y_kx,11));

n3_temp  = size(dat_mb_y_kx,3);
n5_temp  = size(dat_mb_y_kx,5);
n6_temp  = size(dat_mb_y_kx,6);
n7_temp  = size(dat_mb_y_kx,7);
n8_temp  = size(dat_mb_y_kx,8);
n9_temp  = size(dat_mb_y_kx,9);
n10_temp  = size(dat_mb_y_kx,10);
n11_temp  = size(dat_mb_y_kx,11);

if (para.NSeg>1)  && (para.bIntSCorr==1)
    phase_map_temp = phase_map./(abs(phase_map)+1e-10); % Self-navigation
elseif (para.bNav==1)
    phase_map_temp = img_Nav./(abs(img_Nav)+1e-10); phase_map_temp=conj(phase_map_temp); % navigator
else
    phase_map_temp=1;
end

if para.bDisCorr==1
    B0_fmaps_temp = B0_fmaps;
else
    B0_fmaps_temp=1;
end
for idx_RO = 1:ny
    for idx_5 = 1:n5_temp
        for idx_6 = 1:n6_temp
            for idx_7 = 1:n7_temp
                for idx_8 = 1:n8_temp
                    for idx_9 = 1:n9_temp
                        for idx_10 = 1:n10_temp
                            for idx_11 = 1:n11_temp

                                dat_temp = dat_mb_y_kx(idx_RO,:,:,:,idx_5,...
                                    idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);
                                dat_temp = dat_temp(:);

                                temp_sense = maps_sense(idx_RO,:,:,idx_5);
                                temp_sense = repmat(temp_sense,[sum(P==1)*para.NSeg,1,1]);
                                temp_sense = permute(temp_sense,[1,3,2]);
                                temp_sense = reshape(temp_sense,[],nx);

                                if phase_map_temp~=1
                                    phase_temp = phase_map_temp(idx_RO,:,:,idx_5,idx_6,idx_7,idx_8,idx_9,idx_10,idx_11);
                                    phase_temp = repmat(phase_temp,sum(P==1),1,1,size(dat_mb_y_kx,4));
                                    phase_temp = permute(phase_temp,[1,3,4,2]);
                                    phase_temp = reshape(phase_temp,[],size(phase_temp,4));
                                else
                                    phase_temp=phase_map_temp;
                                end

                                if para.bDisCorr==1
                                    B0maps_temp = B0_fmaps_temp(idx_RO,:,idx_5);
                                    B0_scale = (-para.ind_kc+(1:sum(P==1)))./para.deltaTE.*para.ESP;B0_scale = reshape(B0_scale,[],1);
                                    B0maps_temp = exp(1i.*angle(B0maps_temp).*B0_scale);
                                    B0maps_temp = repmat(B0maps_temp,[para.NCha*para.NSeg,1]);
                                else
                                    B0maps_temp=B0_fmaps_temp;
                                end
                                E_sense = E.*temp_sense.*(phase_temp).*B0maps_temp;

                                [x_temp,flag,relres,iter,resvec,lsvec] = lsqr(E_sense,dat_temp,[],para.iterLSQR_combShot);

                                % other solvers possible
                                % x_temp = pinv(E_sense)*dat_temp;

                                img_sense(idx_RO,:,:,idx_5,...
                                    idx_6,idx_7,idx_8,idx_9,idx_10,idx_11) = x_temp;
                            end
                        end
                    end
                end
            end
        end
    end
end


dat_sense = SymIfft(img_sense,[1 2]);
% Eliminate zeros along readout per segments
dat_sense = dat_sense.*reshape(ref_MTF_shots,[],1,size(ref_MTF_shots,2));
% eliminate the zero k-space locations
if para.nPFT~=1
    dat_sense(:,1:(loc_pe-1),:,:,:,:,:,:,:,:,:)=0;
end

ref_MTF(ref_MTF<1e-4)=1;
dat_sense = sum(dat_sense,3);
dat_sense = dat_sense./ref_MTF;
img_sense = SymFft(dat_sense,[1 2]);
dat_epi_comb = SymIfft(img_sense,[1 2]);

%% Apply partial Fourier transform

if para.nPFT ~= 1

    nIter = 8;

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
img_epi_ms_sos=squeeze(abs(img_epi_comb));

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
