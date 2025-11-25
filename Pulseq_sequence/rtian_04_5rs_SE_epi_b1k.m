%% Script description

% This demo script can produce multi-shot readout-segmented sequences,
% adapted from writeEpiDiffusionRS.m in Pulseq'demo folder
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
% 5-shot, readout-segmented, 2-slice, b=0,1000s/mm2, 2D navigator ON

% We suggest a few basic tests:
% 1. produce EPI sequence: set "mode_proj = 0";
% 2. produce separate EPI 1D navigator scan, for N/2 ghost correction: "set mode_proj = 1";

% It is possible to generate EPI with different shot number, diffusion
% weighting, spin-echo/gradient-echo. Please check protocol input. 

% Please ensure to check the acoustic resonance bands at the end of the
% script. Running EPI with improper echo spacing could trigger mechanical oscillation of gradients and
% cause damages.

% It is possible to examine Pulseq sequence visualization and detailed reports by setting mode_vis = 1;

%% User input: please enter for scan protocols

% Trajectory mode
% 0: readout segmented EPI. 
% 1: separate 1D navigator scan with phase encoding disabled.
mode_proj = 0;

% Visualization mode
% 0: Do not exam sequence visualization and reports. 
% 1: visualization and reports. can be slow.
mode_vis = 0;

% Set system limits (Prisma max. 80mT/m@200T/m/s, was using 50T/m/,180T/m/s)
lims = mr.opts('MaxGrad',43,'GradUnit','mT/m',...
    'MaxSlew',140,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89);
lims_refocus_1 = mr.opts('MaxGrad',23,'GradUnit','mT/m',...
    'MaxSlew',110,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89);
lims_refocus_2 = mr.opts('MaxGrad',23,'GradUnit','mT/m',...
    'MaxSlew',110,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89);
lims_diff = mr.opts('MaxGrad',43,'GradUnit','mT/m',...
    'MaxSlew',170,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89);
lims_prephase = mr.opts('MaxGrad',30,'GradUnit','mT/m',...
    'MaxSlew',120,'SlewUnit','T/m/s',...
    'rfRingdownTime', 10e-6, 'rfDeadtime', 100e-6, 'B0', 2.89);
seq=mr.Sequence(lims); % Create a new sequence object

% ------------------------------------------------------
% User input below
% ------------------------------------------------------
gyro_kHzmT = 42.57; % gyromagnetic ratio, gamma/2Pi
fov=220e-3; % Define FOV [mm]
Ny=320; % Define pixel size along phase encoding
Nx=75; % Define pixel size along readout
b2DNav = 1; % 0: no 2D navigator. 1: 2D navigator
Nav_Ny = 40; % 2D navigator: phase encoding pixel size
Nav_Nx = Nx; % 2D navigator: readout pixel size
Nx_op=20; % Define the minimum required overlapped pixel number along readout
Nshot = 5; % Number of readout segments, should be odd for now
bSP = 1; % 1: switching readout polarity between shots. 0: identical readout polarity between shots
bSegShift = 1; % if segment shift is enabled. 0 :disabled, 1: enabled
nPreP = 1+(Nshot-1)/2; % multiplier to adjust the pre phase duration
spoilFactor=1.5; % scale factor of spoiling gradient around the pi-pulse. 0: disabled crusher
bTypeEcho = 1; % 1:spin echo; 2:gradient echo; others: error

if mod(Nshot,2)~=1
    error('even shot number not allowed!')
end
Nx_tot = Nx*Nshot-(Nshot-1)*Nx_op;
if ismember(bTypeEcho,[1,2])==0
    error('echo type is neither spin nor gradient echo!')
end
disp(' ')
disp(['The final pixel size for multi-shot readout segments is: ',num2str(Nx_tot)])
disp(['Multi-shot resolution enhancement factor: ',num2str(round(Nx_tot/Nx,1))])
disp(' ')

thickness = 3e-3; % slice thinckness [m]
dis_sl = 1.5; % slice distance between slice centers (instead of gap as in Siemens), e.g., 2: centers of neighboring slices are 2x slice thickness
Nslices = 2; % slice number
bFactor = 1000; % [s/mm^2]
b3Trace = 0; % If three orthogonal diffusion directions are selected instead of golden angle spiral generation
TE = 85e-3; % TE [s]
TR = 4; % TR [s]
nDummy = 6; % dummy scan number, for spin steady state
nRep = 1; % repetition number, >>1 for fMRI scan
ndd = 2; % diffusion direction number. no diffusion weighted->"0", diffusion weighted direction 64->"64"
nPI = 3; % parallel imaging acceleration factor along phase encoding dimension
vFlip_deg = 58; % RF excitation angle. Here, only allows for GRE scan.

bGzRev = 1; % fat-suppressed gradient modifications, 0: no effects, 1: excitation gradient reverse polarity, 2: excitation pulse duration prolonged
flip_fat = 100; % flip angle of fat saturation [degree]

pe_enable=1;            % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
ro_os=1;                % oversampling factor (in contrast to the product sequence we don't really need it)
readoutTime=4.4e-4;     % this controls the readout bandwidth, standard 3.2e-4
partFourierFactor=6/8;  % partial Fourier factor: 1: full sampling

tBuffTrain = 1e-3; %[s], the buffer between pre-phaser and EPI readout train to reduce DW gradients' eddy currents, e.g. 1ms

% Default RF durations
tRFex=3e-3;
tRFref=3e-3;

% Modified RF durations if fat-supressed gradient modification is used.
if bGzRev==1
    tRFex = 5e-3;
    tRFref = 5e-3;
elseif bGzRev==2
    tRFex = 6.7e-3;
    tRFref= 2.56e-3;
    % tRFex = 13e-3; % 3T condition, but TE will be too longer
end

readoutBW = 1/readoutTime ; % readout bandwidth
disp(['Low-resolution readout bandwidth = ', num2str(readoutBW), ' Hz/Px']) ;

if (bTypeEcho==2) && (ndd>0)
    error('gradient echo is not compatible with diffusion gradient')
end

%% ===========================================================
% Below: sequence generation codes without need of user input
% ============================================================

%% Automatic overwrites, if using 1D navigator scan for N/2 correction
if mode_proj==1
    b2DNav = 0;
    bSP = 0;
    Nslices=1;
    bFactor=0;
    nDummy = 0;
    nRep=1;
    ndd=0;
    pe_enable=0;
    vFlip_deg = 90;
end

%% Construct sequence objects

% Create fat-sat pulse
sat_ppm=-3.3;
sat_freq=sat_ppm*1e-6*lims.B0*lims.gamma;
rf_fs = mr.makeGaussPulse(flip_fat*pi/180,'system',lims,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase
gz_fs = mr.makeTrapezoid('z',lims,'delay',mr.calcDuration(rf_fs),'Area',0.7/1e-4); % spoil up to 0.1mm
gx_fs = mr.makeTrapezoid('x',lims,'delay',mr.calcDuration(rf_fs)+gz_fs.riseTime+0.3e-3,'Area',0.1/1e-4); % spoil up to 0.1mm
gy_fs = mr.makeTrapezoid('y',lims,'delay',mr.calcDuration(rf_fs)+gz_fs.riseTime+gx_fs.riseTime+0.3e-3,'Area',0.3/1e-4); % spoil up to 0.1mm
gz_fs_neg = mr.makeTrapezoid('z',lims,'Area',-0.7/1e-4); % spoil up to 0.1mm
gx_fs_neg = mr.makeTrapezoid('x',lims,'delay',gz_fs_neg.riseTime+0.3e-3,'Area',-0.1/1e-4); % spoil up to 0.1mm
gy_fs_neg = mr.makeTrapezoid('y',lims,'delay',gz_fs_neg.riseTime+gx_fs_neg.riseTime+0.3e-3,'Area',-0.3/1e-4); % spoil up to 0.1mm

% Create RF excitation pulse
if bTypeEcho ==1
    % Create 90 degree slice selection pulse and gradient
    [rf, gz, gzReph] = mr.makeSincPulse(pi/2,'system',lims,'Duration',tRFex,...
        'SliceThickness',thickness,'PhaseOffset',pi/2,'apodization',0.5,'timeBwProduct',4);
else
    % Create alpha degree slice selection pulse and gradient
    [rf, gz, gzReph] = mr.makeSincPulse(vFlip_deg/180*pi,'system',lims,'Duration',tRFex,...
        'SliceThickness',thickness,'PhaseOffset',pi/2,'apodization',0.5,'timeBwProduct',4);
end

% Gradient revsersal for fat suppression
if bGzRev==1
    gz.amplitude = -gz.amplitude;
    gz.area = -gz.area;
    gz.flatArea = -gz.flatArea;
    gzReph.amplitude = -gzReph.amplitude;
    gzReph.area = -gzReph.area;
    gzReph.flatArea = -gzReph.flatArea;
end

% Create 180 degree slice refocusing pulse and gradients
[rf180, gz180] = mr.makeSincPulse(pi,'system',lims_refocus_1,'Duration',tRFref,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'use','refocusing');
if spoilFactor==0
    [~, gzr_t, gzr_a]=mr.makeExtendedTrapezoidArea('z',gz180.amplitude,0,-gzReph.area+0.5*gz180.amplitude*gz180.fallTime,lims_refocus_1);
    gz180n=mr.makeExtendedTrapezoid('z','system',lims_refocus_1,'times',[0 gz180.riseTime gz180.riseTime+gz180.flatTime+gzr_t]+gz180.delay, 'amplitudes', [0 gz180.amplitude gzr_a]);
else
    [~, gzr1_t, gzr1_a]=mr.makeExtendedTrapezoidArea('z',0,gz180.amplitude,spoilFactor*gz180.area,lims_refocus_1);
    [~, gzr2_t, gzr2_a]=mr.makeExtendedTrapezoidArea('z',gz180.amplitude,0,-gzReph.area+spoilFactor*gz180.area,lims_refocus_1);
    if gz180.delay>(gzr1_t(4)-gz180.riseTime)
        gz180.delay=gz180.delay-(gzr1_t(4)-gz180.riseTime);
    else
        rf180.delay=rf180.delay+(gzr1_t(4)-gz180.riseTime)-gz180.delay;
        gz180.delay=0;
    end
    gz180n=mr.makeExtendedTrapezoid('z','system',lims_refocus_1,'times',[gzr1_t gzr1_t(4)+gz180.flatTime+gzr2_t]+gz180.delay, 'amplitudes', [gzr1_a gzr2_a]);
end

% Check if fat-suppressed gradient modifications will work.
if bGzRev==1
    lterm = abs(sat_ppm*lims.B0*42.58/(4/tRFref)); % mT/m
    rterm = 1 - abs(sat_ppm)*lims.B0*42.58/(4/tRFex);
    if lterm>=rterm
        disp('Fat shift supp. condition matched (varied polarity).')
    else
        disp('Fat shift supp. condition NOT matched (varied polarity).')
    end
elseif bGzRev==2
    lterm = abs(gz.amplitude - gz180.amplitude)/42.58e6*1e3; % mT/m
    rterm = 1/abs(2*sat_ppm*lims.B0*42.58)*...
        (4*abs(gz.amplitude/42.58e6*1e3)/tRFref + 4*abs(gz180.amplitude/42.58e6*1e3)/tRFex);
    if lterm>=rterm
        disp('Fat shift supp. condition matched (varied duration).')
    else
        disp('Fat shift supp. condition NOT matched (varied duration).')
    end
end

% Create 2D navigator refocusing pulse
[Nav_rf180, Nav_gz180] = mr.makeSincPulse(pi,'system',lims_refocus_2,'Duration',tRFref,...
    'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'use','refocusing');
if spoilFactor==0
    [~, Nav_gzr_t, Nav_gzr_a]=mr.makeExtendedTrapezoidArea('z',Nav_gz180.amplitude,0,0.5*Nav_gz180.amplitude*Nav_gz180.fallTime,lims_refocus_2);
    Nav_gz180n=mr.makeExtendedTrapezoid('z','system',lims_refocus_2,'times',[0 Nav_gz180.riseTime Nav_gz180.riseTime+Nav_gz180.flatTime+Nav_gzr_t]+Nav_gz180.delay, 'amplitudes', [0 Nav_gz180.amplitude Nav_gzr_a]);
else
    [~, gzr1_t, gzr1_a]=mr.makeExtendedTrapezoidArea('z',0,Nav_gz180.amplitude,spoilFactor*Nav_gz180.area,lims_refocus_2);
    [~, gzr2_t, gzr2_a]=mr.makeExtendedTrapezoidArea('z',Nav_gz180.amplitude,0,spoilFactor*Nav_gz180.area,lims_refocus_2);
    if Nav_gz180.delay>(gzr1_t(4)-Nav_gz180.riseTime)
        Nav_gz180.delay=Nav_gz180.delay-(gzr1_t(4)-Nav_gz180.riseTime);
    else
        Nav_rf180.delay=Nav_rf180.delay+(gzr1_t(4)-Nav_gz180.riseTime)-Nav_gz180.delay;
        Nav_gz180.delay=0;
    end
    Nav_gz180n=mr.makeExtendedTrapezoid('z','system',lims_refocus_2,'times',[gzr1_t gzr1_t(4)+Nav_gz180.flatTime+gzr2_t]+Nav_gz180.delay, 'amplitudes', [gzr1_a gzr2_a]);
end

% define the output trigger to play out with every slice excitatuion
trig=mr.makeDigitalOutputPulse('osc0','duration', 100e-6); % possible channels: 'osc0','osc1','ext1'

% Define other gradients and ADC events
deltak=1/fov;
kWidth = Nx*deltak;
k_Nav_Width = Nav_Nx*deltak;

% Phase blip in shortest possible time
blip_dur = ceil(2*sqrt(nPI*deltak/lims.maxSlew)/10e-6/2)*10e-6*2; % we round-up the duration to 2x the gradient raster time
% the split code below fails if this really makes a trpezoid instead of a triangle...
gy = mr.makeTrapezoid('y',lims,'Area',-deltak*nPI,'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center
%gy = mr.makeTrapezoid('y',lims,'amplitude',deltak/blip_dur*2,'riseTime',blip_dur/2, 'flatTime', 0);

% readout gradient is a truncated trapezoid with dead times at the beginnig
% and at the end each equal to a half of blip_dur
% the area between the blips should be defined by kWidth
% we do a two-step calculation: we first increase the area assuming maximum
% slewrate and then scale down the amlitude to fix the area
extra_area=blip_dur/2*blip_dur/2*lims.maxSlew; % check unit!; % rtian assume the shortest triangle at sides of rectangular gradient pulse
gx = mr.makeTrapezoid('x',lims,'Area',kWidth+extra_area,'duration',readoutTime+blip_dur);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx.amplitude=gx.amplitude/actual_area*kWidth; % rtian: re-scale the gradient amplitude
gx.area = gx.amplitude*(gx.flatTime + gx.riseTime/2 + gx.fallTime/2);
gx.flatArea = gx.amplitude*gx.flatTime;
ESP = 1e3 * mr.calcDuration(gx) ; % echo spacing, ms
disp(['echo spacing = ', num2str(ESP), ' ms']) ;

% The readout gradient replicated for navigator
gnav_x = mr.makeTrapezoid('x',lims,'Area',k_Nav_Width+extra_area,'duration',readoutTime*Nav_Nx/Nx+blip_dur);
actual_nav_area=gnav_x.area-gnav_x.amplitude/gnav_x.riseTime*blip_dur/2*blip_dur/2/2-gnav_x.amplitude/gnav_x.fallTime*blip_dur/2*blip_dur/2/2;
gnav_x.amplitude=gnav_x.amplitude/actual_nav_area*k_Nav_Width; % rtian: re-scale the gradient amplitude
gnav_x.area = gnav_x.amplitude*(gnav_x.flatTime + gnav_x.riseTime/2 + gnav_x.fallTime/2);
gnav_x.flatArea = gnav_x.amplitude*gnav_x.flatTime;

% Calculate ADC
% we use ramp sampling, so we have to calculate the dwell time and the
% number of samples, which are will be qite different from Nx and
% readoutTime/Nx, respectively.
adcDwellNyquist=deltak/gx.amplitude/ro_os;
% round-down dwell time to 100 ns
adcDwell=floor(adcDwellNyquist*1e7)*1e-7;
adcSamples=floor(readoutTime/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
% MZ: no idea, whether ceil,round or floor is better for the adcSamples...
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',blip_dur/2);
% realign the ADC with respect to the gradient
time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % I've been told that Siemens samples in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)*1e6)*1e-6; % we adjust the delay to align the trajectory with the gradient. We have to aligh the delay to 1us
% this rounding actually makes the sampling points on odd and even readouts
% to appear misalligned. However, on the real hardware this misalignment is
% much stronger anyways due to the grdient delays

% Calculate ADC for navigator
Nav_adcDwellNyquist=deltak/gnav_x.amplitude/ro_os;
Nav_adcDwell=floor(Nav_adcDwellNyquist*1e7)*1e-7;
Nav_adcSamples=floor(readoutTime*Nav_Nx/Nx/Nav_adcDwell/4)*4;
Nav_adc = mr.makeAdc(Nav_adcSamples,'Dwell',Nav_adcDwell,'Delay',blip_dur/2);
Nav_time_to_center=Nav_adc.dwell*((Nav_adcSamples-1)/2+0.5);
Nav_adc.delay=round((gnav_x.riseTime+gnav_x.flatTime/2-Nav_time_to_center)*1e6)*1e-6;

% FOV positioning requires alignment to grad. raster... -> TODO

% Split the blip into two halves and produce a combined synthetic gradient
gy_parts = mr.splitGradientAt(gy, blip_dur/2, lims);
[gy_blipup, gy_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gx);
gy_blipdownup=mr.addGradients({gy_blipdown, gy_blipup}, lims);

% Navigator: split the blip
[gnav_y_blipup, gnav_y_blipdown,~]=mr.align('right',gy_parts(1),'left',gy_parts(2),gnav_x);
gnav_y_blipdownup=mr.addGradients({gnav_y_blipdown, gnav_y_blipup}, lims);

% pe_enable support
gy_blipup.waveform=gy_blipup.waveform*pe_enable;
gy_blipdown.waveform=gy_blipdown.waveform*pe_enable;
gy_blipdownup.waveform=gy_blipdownup.waveform*pe_enable;
gnav_y_blipup.waveform = gnav_y_blipup.waveform*pe_enable;
gnav_y_blipdown.waveform = gnav_y_blipdown.waveform*pe_enable;
gnav_y_blipdownup.waveform = gnav_y_blipdownup.waveform*pe_enable;

% phase encoding and partial Fourier
Ny_pre=round((1/nPI)*partFourierFactor*Ny/2-1); % PE steps prior to ky=0, excluding the central line
Ny_post=round((1/nPI)*Ny/2); % PE lines after the k-space center including the central line, the Maxim's +1 is deleted to ensure 0 @ center
Ny_meas=Ny_pre+Ny_post;
disp(['Phase encoded lines per shot: ',num2str(Ny_meas)])
disp(['Acq. window about: ',num2str((Ny_meas+1)*ESP),' ms'])
disp(['The central k-space line is at #',num2str(Ny_pre+1), ' phase encoded lines'])

% phase encoding for navigator
Ny_nav_pre=round((1/nPI)*Nav_Ny/2-1); % PE steps prior to ky=0, excluding the central line
Ny_nav_post=round((1/nPI)*Nav_Ny/2); % PE lines after the k-space center including the central line
Ny_nav_meas = Ny_nav_pre+Ny_nav_post;
disp(['The central k-space line for 2D nav is at #',num2str(Ny_nav_pre+1), ' phase encoded lines'])
disp(' ')

% Pre-phasing gradients
nsb = (Nshot-1)/2;
gxPre = mr.makeTrapezoid('x',lims_prephase,'Area',-gx.area*(1/2+nsb),'duration',ceil(nPreP*readoutTime./lims.gradRasterTime).*lims.gradRasterTime);
gyPre = mr.makeTrapezoid('y',lims_prephase,'Area',Ny_pre*deltak*nPI,'duration',ceil(nPreP*readoutTime./lims.gradRasterTime).*lims.gradRasterTime);

[gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
% relax the PE prepahser to reduce stimulation
gyPre = mr.makeTrapezoid('y',lims_prephase,'Area',gyPre.area,'Duration',mr.calcDuration(gxPre,gyPre));
gyPre.amplitude=gyPre.amplitude*pe_enable;

% Navigator: pre-phase gradients
Ny_diff = Ny_post-Ny_nav_pre-1;
Nav_gxPre = mr.makeTrapezoid('x',lims_prephase,'Area',-gnav_x.area/2,'duration',ceil(nPreP*readoutTime*Nav_Nx/Nx./lims.gradRasterTime).*lims.gradRasterTime);
Nav_gyPre = mr.makeTrapezoid('y',lims_prephase,'Area',-Ny_diff*deltak*nPI,'duration',ceil(nPreP*readoutTime*Nav_Nx/Nx./lims.gradRasterTime).*lims.gradRasterTime);

[Nav_gxPre,Nav_gyPre]=mr.align('right',Nav_gxPre,'left',Nav_gyPre);
% relax the PE prepahser to reduce stimulation
Nav_gyPre = mr.makeTrapezoid('y',lims_prephase,'Area',Nav_gyPre.area,'Duration',mr.calcDuration(Nav_gxPre,Nav_gyPre));
Nav_gyPre.amplitude=Nav_gyPre.amplitude*pe_enable;

% Calculate delay times
durationToCenter = (Ny_pre+0.5)*mr.calcDuration(gx); % rtian: start of EPI readout to sampling the k-space ceter
rfCenterInclDelay=rf.delay + mr.calcRfCenter(rf); % rtian: half of RF block
    
if bTypeEcho ==1 % spin echo
    rf180centerInclDelay=rf180.delay + mr.calcRfCenter(rf180);% rtian: half of RF block
    delayTE1=ceil((TE/2 - mr.calcDuration(rf,gz) + rfCenterInclDelay - rf180centerInclDelay)/lims.gradRasterTime)*lims.gradRasterTime;
    delayTE2tmp=ceil((TE/2 - mr.calcDuration(rf180,gz180n) + rf180centerInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
    assert(delayTE1>=0);
else % gradient echo, take delayTR2, discard delayTE1
    delayTE2tmp=ceil((TE - mr.calcDuration(rf,gz) -  mr.calcDuration(gzReph) + rfCenterInclDelay - durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
end
%delayTE2=delayTE2tmp+mr.calcDuration(rf180,gz180n);
gxPre.delay=0;
gyPre.delay=0;
delayTE2=delayTE2tmp-mr.calcDuration(gxPre,gyPre)-tBuffTrain;
[gxPre,gyPre]=mr.align('right',gxPre,'left',gyPre);
assert(delayTE2>=0);

% Navigator: calculate delay time based on TE 2
Nav_durationToCenter = (Ny_nav_pre+0.5)*mr.calcDuration(gnav_x);
CenterToEnd = (Ny_post-0.5)*mr.calcDuration(gx);
Nav_rf180centerInclDelay = Nav_rf180.delay + mr.calcRfCenter(Nav_rf180);% rtian: half of RF block;
delayTE3 = ceil((CenterToEnd + mr.calcDuration(gxPre) + Nav_rf180centerInclDelay - (mr.calcDuration(Nav_rf180,Nav_gz180n)-Nav_rf180centerInclDelay) - mr.calcDuration(Nav_gyPre) - tBuffTrain - Nav_durationToCenter)/lims.gradRasterTime)*lims.gradRasterTime;
TE2_nav = TE+CenterToEnd+mr.calcDuration(gxPre)+mr.calcDuration(Nav_rf180,Nav_gz180n)+ mr.calcDuration(Nav_gyPre)+ delayTE3+ tBuffTrain + Nav_durationToCenter;

if b2DNav==1
    assert(delayTE3>=0);
    disp(['TE 2D navigator = ', num2str(TE2_nav*1e3), ' ms']) ;
end

% Calculate diffusion gradients
if bTypeEcho == 1
    % diffusion weithting calculation
    % delayTE2 is our window for small_delta
    % delayTE1+delayTE2-delayTE2 is our big delta
    % we anticipate that we will use the maximum gradient amplitude, so we need
    % to shorten delayTE2 by gmax/max_sr to accommodate the ramp down
    gr=ceil(lims_diff.maxGrad/lims_diff.maxSlew/lims_diff.gradRasterTime)*lims_diff.gradRasterTime;% rtian
    small_delta=delayTE2-gr;%rtian
    big_delta=delayTE1+mr.calcDuration(rf180,gz180n);
    % we define bFactCalc function below to eventually calculate time-optimal
    % gradients. for now we just abuse it with g=1 to give us the coefficient
    epi = gr; % rtian: for trapezoid
    g =sqrt(bFactor*1e6/bFactCalc(1,small_delta,big_delta,epi)); % for now it looks too large!
    % rtian: seems like g is in [Hz/m]
    
    gDiff_x=mr.makeTrapezoid('x','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims_diff);
    gDiff_y=mr.makeTrapezoid('y','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims_diff);
    gDiff_z=mr.makeTrapezoid('z','amplitude',g,'riseTime',gr,'flatTime',small_delta-gr,'system',lims_diff);
    
    assert(mr.calcDuration(gDiff_x)<=delayTE1);
    assert(mr.calcDuration(gDiff_x)<=delayTE2);
    assert(mr.calcDuration(gDiff_y)<=delayTE1);
    assert(mr.calcDuration(gDiff_y)<=delayTE2);
    assert(mr.calcDuration(gDiff_z)<=delayTE1);
    assert(mr.calcDuration(gDiff_z)<=delayTE2);
    
end

% Spoiler gradient, not used now
factor_sp = 6;
gy_spoiler = mr.makeTrapezoid('y',lims_prephase,'Area',-factor_sp*Ny_post*deltak*nPI,'duration',readoutTime/Nx*Ny_post*nPI*factor_sp);
gz_spoiler = mr.makeTrapezoid('z',lims_prephase,'Area',factor_sp*Ny_post*deltak*nPI,'duration',readoutTime/Nx*Ny_post*nPI*factor_sp);
[gy_spoiler,gz_spoiler]=mr.align('right',gy_spoiler,'right',gz_spoiler);
gy_spoiler.flatTime = round(gy_spoiler.flatTime*1e6,-1)/1e6;
gz_spoiler.flatTime = round(gz_spoiler.flatTime*1e6,-1)/1e6;

% For even/odd echoes collection
if mode_proj==1
    Ny_meas = 50;
    Nshot =1;
end

disp(['readout gradient amplitude = ', num2str(gx.amplitude/42.58e6*1e3), ' mT/m']) ;
disp(['readout slew rate = ', num2str(gx.amplitude/42.58e6/gx.riseTime), ' T/m/s']) ;
disp(' ')

%% Dummy scan

% Define dummy scans
for idx_ds = 1:nDummy
    ms=1; % always run the 1st segment
    for s=[(1:2:Nslices),(2:2:Nslices)]
        
        % RF excitation
        seq.addBlock(rf_fs,gz_fs,gx_fs,gy_fs); %rtian fat sat
        rf.freqOffset=gz.amplitude*thickness*(-dis_sl)*(s-1-(Nslices-1)/2);
        rf.phaseOffset=pi/2-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
        if bTypeEcho == 1
            rf180.freqOffset=gz180.amplitude*thickness*(-dis_sl)*(s-1-(Nslices-1)/2);
            rf180.phaseOffset=-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase
        end
        seq.addBlock(rf,gz,trig);
        
        if bTypeEcho == 1
            % Diffusion gradient
            seq.addBlock(mr.makeDelay(delayTE1));
            
            % 180 refocusing pulse
            seq.addBlock(rf180,gz180n);
            
            % Diffusion gradient
            seq.addBlock(mr.makeDelay(delayTE2));
        else % only delay
            seq.addBlock(gzReph);
            seq.addBlock(mr.makeDelay(delayTE2));
        end
        
        % Calculate prephase readout gradient scale factor
        ind_band = -(ms-1-(Nshot-1)/2);
        
        if (bSP==1) && (mod(ind_band,2)==1)
            gx.amplitude = -gx.amplitude;
            if bSegShift==1
                area_temp = -gx.area*(1/2+ind_band-Nx_op/Nx*ind_band-1);
            else
                area_temp = -gx.area*(-1/2);
            end
            factor = area_temp/(-gx.area*(1/2+nsb));
        else
            if bSegShift==1
                area_temp = -gx.area*(1/2+ind_band-Nx_op/Nx*ind_band);
            else
                area_temp = -gx.area*(1/2);
            end
            factor = area_temp/(-gx.area*(1/2+nsb));
        end
        
        % Create and run the scaled readout prephaser
        gxPre_temp = gxPre;
        [gxPre_temp,t_red] = scaleGrad(lims,gxPre_temp,factor);
        if t_red<(-1e-12)
            error('The gradient cannot be prolonged for timing!')
        end
        
        [gxPre_temp,gyPre]=mr.align('right',gxPre_temp,'right',gyPre);
        seq.addBlock(gyPre,gxPre_temp);
        
        % buffer to avoid eddy current
        if tBuffTrain>0
            seq.addBlock(mr.makeDelay(tBuffTrain));
        end
        
        % imaging echoes
        for i=1:Ny_meas
            if i==1
                seq.addBlock(gx,gy_blipup); % Read the first line of k-space with a single half-blip at the end
            elseif i==Ny_meas
                seq.addBlock(gx,gy_blipdown); % Read the last line of k-space with a single half-blip at the beginning
            else
                seq.addBlock(gx,gy_blipdownup); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
            end
            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
        end

        % If odd PE lines, flip the gx to default
        if mod(Ny_meas,2)==1
            gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
        end

        % If the band readout polarity is flipped, flip the gx
        if (bSP==1) && (mod(ind_band,2)==1)
            gx.amplitude = -gx.amplitude;
        end

        % 2D navigator echoes
        if b2DNav == 1
            gNav_xPre_temp = gxPre;
            oddness_phase = mod(Ny_meas,2);
            if (bSP==1) && (mod(ind_band,2)==1)
                sign_control = -(0.5 + oddness_phase);
            else
                sign_control = (-1)^(oddness_phase+1)*0.5;
            end

            % Scale the readout prephase for 2D navigator
            [gNav_xPre_temp,t_red] = scaleGrad(lims,gNav_xPre_temp,-factor+sign_control/(1/2+nsb));

            % Run the readout prephase for 2D navigator
            seq.addBlock(gNav_xPre_temp);
            if t_red>0
                seq.addBlock(mr.makeDelay(t_red));
            elseif t_red<(-1e-12)
                error('The gradient cannot be prolonged for timing!')
            end

            % 180 pulse slice selective
            Nav_rf180.freqOffset=Nav_gz180.amplitude*thickness*(-dis_sl)*(s-1-(Nslices-1)/2);
            Nav_rf180.phaseOffset=-2*pi*Nav_rf180.freqOffset*mr.calcRfCenter(Nav_rf180); % compensate for the slice-offset induced phase

            % Run the 180 pulse
            seq.addBlock(Nav_rf180,Nav_gz180n);

            % Run the phase prephase for 2D navigator
            seq.addBlock(Nav_gyPre);

            % Run the delay for 2D navigator. If buffer large,
            % this becomes small, as in calculation
            seq.addBlock(mr.makeDelay(delayTE3));

            % buffer to avoid eddy current
            if tBuffTrain>0
                seq.addBlock(mr.makeDelay(tBuffTrain));
            end

            % 2D Navigator
            for i=1:Ny_nav_meas
                if i==1
                    seq.addBlock(gnav_x,gnav_y_blipup); % Read the first line of k-space with a single half-blip at the end
                elseif i==Ny_nav_meas
                    seq.addBlock(gnav_x,gnav_y_blipdown); % Read the last line of k-space with a single half-blip at the beginning
                else
                    seq.addBlock(gnav_x,gnav_y_blipdownup); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                end
                gnav_x.amplitude = -gnav_x.amplitude;   % Reverse polarity of read gradient
            end
            if mod(Ny_nav_meas,2)==1 % If odd PE lines, flip the gx to default
                gnav_x.amplitude = -gnav_x.amplitude;   % Reverse polarity of read gradient
            end

        end

        % Add spoiler
        seq.addBlock(gz_fs_neg,gy_fs_neg,gx_fs_neg); %rtian fat sat

        % Fill the rest ilde time given the fixed TR
        if (s==1) && (idx_ds==1)
            ts_block = sum(seq.blockDurations);
            ts_idle = TR/Nslices-ts_block;
            if ts_idle<0
                error('too many slices')
            else
                disp(['scan time for one acquisition = ', num2str(ts_block*1e3), ' ms']) ;
            end
        end
        
        seq.addBlock(mr.makeDelay(round(ts_idle/lims.gradRasterTime)*lims.gradRasterTime));
    end
end



%% MR Acquisition

% Create diffusion direction vector
if ndd>=1 && (b3Trace==0)
    dir_diff = rtian_GenDiffusionDirection(ndd);
elseif b3Trace==1
    dir_diff = [0,0,1;0,1,0;1,0,0];
end

% Define sequence blocks
seq.addBlock(mr.makeLabel('SET','REP', 0)); % remove according to Maxim,
% so that the number of run in exam card can be effective.
for idx_avg = 1:nRep
    
    seq.addBlock(mr.makeLabel('SET','PHS',0));
    for idx_dd = 1:(ndd+1)
        
        if bTypeEcho == 1
            % Set the diffusion gradient direction
            gDiff_x_temp = gDiff_x;
            gDiff_y_temp = gDiff_y;
            gDiff_z_temp = gDiff_z;
            
            if idx_dd==1 % b=0 image
                gDiff_x_temp.amplitude = 0;
                gDiff_y_temp.amplitude = 0;
                gDiff_z_temp.amplitude = 0;
            else % diffusion weighted image
                gDiff_x_temp.amplitude = gDiff_x_temp.amplitude*dir_diff(idx_dd-1,1);
                gDiff_y_temp.amplitude = gDiff_y_temp.amplitude*dir_diff(idx_dd-1,2);
                gDiff_z_temp.amplitude = gDiff_z_temp.amplitude*dir_diff(idx_dd-1,3);
            end
        end
        
        seq.addBlock(mr.makeLabel('SET','SEG',0)) ;
        for ms = 1:Nshot
            
            seq.addBlock(mr.makeLabel('SET','SET',0)) ; % Change from SLC to SET to avoid Pulseq FOV shift error
            for s=[(1:2:Nslices),(2:2:Nslices)]
                
                % RF excitation
                seq.addBlock(rf_fs,gz_fs,gx_fs,gy_fs); %rtian fat sat
                rf.freqOffset=gz.amplitude*thickness*(-dis_sl)*(s-1-(Nslices-1)/2);
                rf.phaseOffset=pi/2-2*pi*rf.freqOffset*mr.calcRfCenter(rf); % compensate for the slice-offset induced phase
                
                if bTypeEcho == 1
                    rf180.freqOffset=gz180.amplitude*thickness*(-dis_sl)*(s-1-(Nslices-1)/2);
                    rf180.phaseOffset=-2*pi*rf180.freqOffset*mr.calcRfCenter(rf180); % compensate for the slice-offset induced phase
                end
                seq.addBlock(rf,gz,trig);
                
                if bTypeEcho == 1
                    % Diffusion gradient
                    seq.addBlock(mr.makeDelay(delayTE1),gDiff_x_temp,gDiff_y_temp,gDiff_z_temp);
                    
                    % 180 refocusing pulse
                    seq.addBlock(rf180,gz180n);
                    
                    % Diffusion gradient
                    seq.addBlock(mr.makeDelay(delayTE2),gDiff_x_temp,gDiff_y_temp,gDiff_z_temp);
                else % only delay
                    seq.addBlock(gzReph);
                    seq.addBlock(mr.makeDelay(delayTE2));
                end
                
                % Calculate prephase readout gradient scale factor
                ind_band = -(ms-1-(Nshot-1)/2);
                
                if (bSP==1) && (mod(ind_band,2)==1)
                    gx.amplitude = -gx.amplitude;
                    if bSegShift==1
                        area_temp = -gx.area*(1/2+ind_band-Nx_op/Nx*ind_band-1);
                    else
                        area_temp = -gx.area*(-1/2);
                    end
                    factor = area_temp/(-gx.area*(1/2+nsb));
                else
                    if bSegShift==1
                        area_temp = -gx.area*(1/2+ind_band-Nx_op/Nx*ind_band);
                    else
                        area_temp = -gx.area*(1/2);
                    end
                    factor = area_temp/(-gx.area*(1/2+nsb));
                end
                
                % Create and run the scaled readout prephaser
                gxPre_temp = gxPre;
                
                % imaging acquisition
                [gxPre_temp,t_red] = scaleGrad(lims,gxPre_temp,factor);
                if t_red<(-1e-12)
                    error('The gradient cannot be prolonged for timing!')
                end
                
                [gxPre_temp,gyPre]=mr.align('right',gxPre_temp,'right',gyPre);
                seq.addBlock(gyPre,gxPre_temp);
                
                if (idx_avg==1) && (idx_dd==1) && (s==1)
                    disp(['prephase x gradient amplitude = ', num2str(gxPre_temp.amplitude/42.58e6*1e3), ' mT/m']) ;
                    disp(['prephase x slew rate = ', num2str(gxPre_temp.amplitude/42.58e6/gxPre_temp.riseTime), ' T/m/s']) ;
                end

                % buffer to avoid eddy current
                if tBuffTrain>0
                    seq.addBlock(mr.makeDelay(tBuffTrain));
                end
                
                % Start acquisition index counting
                seq.addBlock(mr.makeLabel('SET','NAV', 0),mr.makeLabel('SET','LIN', 0));
                
                % imaging echoes
                for i=1:Ny_meas
                    mr.makeLabel('SET','REV', gx.amplitude<0);
                    if i==1
                        seq.addBlock(gx,gy_blipup,adc); % Read the first line of k-space with a single half-blip at the end
                    elseif i==Ny_meas
                        seq.addBlock(gx,gy_blipdown,adc); % Read the last line of k-space with a single half-blip at the beginning
                    else
                        seq.addBlock(gx,gy_blipdownup,adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                    end
                    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                    seq.addBlock(mr.makeLabel('INC','LIN', 1));
                end
                
                % If odd PE lines, flip the gx to default
                if mod(Ny_meas,2)==1
                    gx.amplitude = -gx.amplitude;   % Reverse polarity of read gradient
                end

                % If the band readout polarity is flipped, flip the gx
                if (bSP==1) && (mod(ind_band,2)==1)
                    gx.amplitude = -gx.amplitude;
                end

                % 2D navigator echoes
                if b2DNav == 1
                    gNav_xPre_temp = gxPre;
                    oddness_phase = mod(Ny_meas,2);
                    if (bSP==1) && (mod(ind_band,2)==1)
                        sign_control = -(0.5 + oddness_phase);
                    else
                        sign_control = (-1)^(oddness_phase+1)*0.5;
                    end

                    % Scale the readout prephase for 2D navigator
                    [gNav_xPre_temp,t_red] = scaleGrad(lims,gNav_xPre_temp,-factor+sign_control/(1/2+nsb));

                    % Run the readout prephase for 2D navigator
                    seq.addBlock(gNav_xPre_temp);
                    if t_red>0
                        seq.addBlock(mr.makeDelay(t_red));
                    elseif t_red<(-1e-12)
                        error('The gradient cannot be prolonged for timing!')
                    end

                    % 180 pulse slice selective
                    Nav_rf180.freqOffset=Nav_gz180.amplitude*thickness*(-dis_sl)*(s-1-(Nslices-1)/2);
                    Nav_rf180.phaseOffset=-2*pi*Nav_rf180.freqOffset*mr.calcRfCenter(Nav_rf180); % compensate for the slice-offset induced phase

                    % Run the 180 pulse
                    seq.addBlock(Nav_rf180,Nav_gz180n);

                    % Run the phase prephase for 2D navigator
                    seq.addBlock(Nav_gyPre);

                    % Run the delay for 2D navigator. If buffer large,
                    % this becomes small, as in calculation
                    seq.addBlock(mr.makeDelay(delayTE3));

                    % buffer to avoid eddy current
                    if tBuffTrain>0
                        seq.addBlock(mr.makeDelay(tBuffTrain));
                    end

                    % Start index counting for 2D navigator
                    seq.addBlock(mr.makeLabel('SET','NAV', 1),mr.makeLabel('SET','LIN', 0));

                    % 2D Navigator
                    for i=1:Ny_nav_meas
                        if i==1
                            seq.addBlock(gnav_x,gnav_y_blipup,Nav_adc); % Read the first line of k-space with a single half-blip at the end
                        elseif i==Ny_nav_meas
                            seq.addBlock(gnav_x,gnav_y_blipdown,Nav_adc); % Read the last line of k-space with a single half-blip at the beginning
                        else
                            seq.addBlock(gnav_x,gnav_y_blipdownup,Nav_adc); % Read an intermediate line of k-space with a half-blip at the beginning and a half-blip at the end
                        end
                        gnav_x.amplitude = -gnav_x.amplitude;   % Reverse polarity of read gradient
                        seq.addBlock(mr.makeLabel('INC','LIN', 1));
                    end
                    if mod(Ny_nav_meas,2)==1 % If odd PE lines, flip the gx to default
                        gnav_x.amplitude = -gnav_x.amplitude;   % Reverse polarity of read gradient
                    end

                end

                % Add spoiler
                seq.addBlock(gz_fs_neg,gy_fs_neg,gx_fs_neg); %rtian fat sat
                
                % Fill the rest ilde time given the fixed TR
                if ((nDummy==0) && (s==1) && (ms ==1) && (idx_dd==1))
                    ts_block = sum(seq.blockDurations);
                    ts_idle = TR/Nslices-ts_block;
                    if ts_idle<0
                        error('too many slices')
                    else
                        disp(['scan time for one acquisition = ', num2str(ts_block*1e3), ' ms']) ;
                    end
                end
                seq.addBlock(mr.makeDelay(round(ts_idle/lims.gradRasterTime)*lims.gradRasterTime));
                
                seq.addBlock(mr.makeLabel('INC','SET', 1)) ; % change from SLC to SET to avoid Pulseq FOV shift error
            end
            seq.addBlock(mr.makeLabel('INC','SEG', 1)) ;
        end
        seq.addBlock(mr.makeLabel('INC','PHS', 1)) ;
    end
    seq.addBlock(mr.makeLabel('INC','REP', 1)) ;
end


%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% do some visualizations

if mode_vis~=0 % could be too slow for large sequence file
    
    seq.plot();             % Plot sequence waveforms
    
    % trajectory calculation
    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
    
    % plot k-spaces
    figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
    hold on; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
    
    figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
    axis('equal'); % enforce aspect ratio for the correct trajectory display
    hold on;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
    
end

%% prepare the sequence output for the scanner
seq.setDefinition('FOV', [fov fov thickness*Nslices]);
seq.setDefinition('Name', 'epi-diff');

seq.write('rtian_epi.seq');

% seq.install('siemens');

% seq.sound(); % simulate the seq's tone

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits
if mode_vis~=0
    rep = seq.testReport;
    fprintf([rep{:}]);
end

gradSpectrum_pulseq_fct_MPITUB;

toc(tst_tot);

%% Diffusion coefficient calculation

function b=bFactCalc(g, delta, DELTA,epi)
% see DAVY SINNAEVE Concepts in Magnetic Resonance Part A, Vol. 40A(2) 3965 (2012) DOI 10.1002/cmr.a
% b = gamma^2  g^2 delta^2 sigma^2 (DELTA + 2 (kappa - lambda) delta)
% in pulseq we don't need gamma as our gradinets are Hz/m
% however, we do need 2pi as diffusion equations are all based on phase
% for rect gradients: sigma=1 lambda=1/2 kappa=1/3
% for trapezoid gradients: TODO
sigma=1;
%lambda=1/2;
%kappa=1/3;
kappa_minus_lambda=1/3-1/2;
% b= (2*pi * g * delta * sigma)^2 * (DELTA + 2*kappa_minus_lambda*delta);

% rtian, p279, Handbook of MRI pusle sequence
g = g*(2*pi);
% delta = delta*1e3;
% DELTA = DELTA*1e3;
% epi = epi*1e3;
b = g^2*delta^2*(DELTA-delta/3) + g^2*epi^3/30 - g^2*delta*epi^2/6; % rtian, trapezoid

end

%% Gradient scaling mostly by duration

function [obj_grad,time_reduced] = scaleGrad(lims,obj_grad,scale)

if abs(scale)>1
    obj_grad.amplitude = obj_grad.amplitude*scale;
    obj_grad.area = obj_grad.area*scale;
    obj_grad.flatArea = obj_grad.flatArea*scale;
    time_reduced = 0;
elseif scale~=0
    area_temp = obj_grad.area*scale;
    flatTime_ori = obj_grad.flatTime;
    if obj_grad.amplitude==0
        error('Gradient amplitude zero! numerical error due to dividing by 0')
    end
    delta_flat = obj_grad.area*(1-abs(scale))/obj_grad.amplitude;
    if obj_grad.flatTime>delta_flat % a reduced trapezoidal
        obj_grad.flatTime = round((obj_grad.flatTime-delta_flat)/lims.gradRasterTime)*lims.gradRasterTime;
    else % a triangle
        obj_grad.flatTime = 0;
    end
    obj_grad.amplitude = obj_grad.amplitude*sign(scale);
    area_actual = (obj_grad.flatTime*2+obj_grad.riseTime+obj_grad.fallTime)*obj_grad.amplitude/2;
    obj_grad.amplitude = obj_grad.amplitude/area_actual*area_temp;
    obj_grad.area=area_temp;
    obj_grad.flatArea = obj_grad.amplitude*obj_grad.flatTime;
    time_reduced = flatTime_ori-obj_grad.flatTime;
else
    obj_grad.amplitude = 0;
    obj_grad.flatArea = 0;
    obj_grad.area = 0;
    time_reduced = 0;
end
end
%% This function calculates diffusion directions based on direction number

function [directions] = rtian_GenDiffusionDirection(num_points)
% cite "A better way to construct the sunflower head"
% Fibonacci sphere point picking
theta = pi * (3. - sqrt(5.));  % golden angle in radians

z = linspace(1,- 1,num_points); % another expression for cos(theta)
radius = sqrt(1 - z.^2); % sin(theta)
phi = theta * (0:num_points-1);

% Spherical to Cartesian coordinates conversion
y = radius .* sin(phi);
x = radius .* cos(phi);

directions = [x', y', z'];

figure
plot3(directions(:,1),directions(:,2),directions(:,3),'o')
title('diffusion weighted directions')

end

