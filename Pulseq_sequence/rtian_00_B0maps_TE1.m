%% Script description

% This demo script can produce GRE sequence for B0 mapping,
% adapted from writeGradientEcho.m in Pulseq'demo folder
% tested with a Siemens 3T Prisma scanner at the Max Planck Institute for Biological Cybernetics in Tuebingen, Germany.

% For questions, please contact: Rui Tian, rui.tian@tuebingen.mpg.de/ruitianwater@outlook.com

addpath(genpath('./'))

%% Sequence

% set system limits
sys = mr.opts('MaxGrad', 40, 'GradUnit', 'mT/m', ...
    'MaxSlew', 120, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq=mr.Sequence(sys);           % Create a new sequence object
fov=220e-3; Nx=300; Ny=300;     % Define FOV and resolution
alpha=10;                       % flip angle
sliceThickness=3e-3;            % slice
TR=20e-3;                       % repetition time TR
TE=8e-3;                        % echo time TE
nDummy = 10;

% slice position
Nslices = 2;
dis_sl = 1.5;

% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
roDuration=ceil(3.2e-3./100.*Nx./seq.gradRasterTime).*seq.gradRasterTime;              % ADC duration

% Create fat-sat pulse
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
B0=2.89; % 1.5 2.89 3.0
sat_ppm=-3.3;
sat_freq=sat_ppm*1e-6*B0*sys.gamma;
rf_fs = mr.makeGaussPulse(100*pi/180,'system',sys,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create alpha-degree slice selection pulse and gradient
% [rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
%     'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);
% apodization ->0.5, compatible with EPI
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
    'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',1e-3,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',1e-3,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
gyPre = mr.makeTrapezoid('y','Area',max(abs(phaseAreas)),'Duration',mr.calcDuration(gxPre),'system',sys);
peScales=phaseAreas/gyPre.area;


% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTE>=0));
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

% Loop over phase encodes and define sequence blocks
for s=1:Nslices
    rf_phase=0;
    rf_inc=0;
    for i=nDummy:-1:1
        for c=1:length(TE)
            seq.addBlock(rf_fs,gz_fs); % fat-sat

            rf.freqOffset=gz.amplitude*sliceThickness*(-dis_sl)*(s-1-(Nslices-1)/2);

            rf.phaseOffset=rf_phase/180*pi;
            adc.phaseOffset=rf_phase/180*pi;
            rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase=mod(rf_phase+rf_inc, 360.0);
            %
            seq.addBlock(rf,gz);
            seq.addBlock(gxPre,mr.scaleGrad(gyPre,peScales(i)),gzReph);
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gx);
            %             gyPre.amplitude=-gyPre.amplitude;
            seq.addBlock(mr.makeDelay(delayTR(c)),gxSpoil,mr.scaleGrad(gyPre,-peScales(i)),gzSpoil)
        end
    end

    for i=Ny:-1:1
        for c=1:length(TE)
            seq.addBlock(rf_fs,gz_fs); % fat-sat

            rf.freqOffset=gz.amplitude*sliceThickness*(-dis_sl)*(s-1-(Nslices-1)/2);

            rf.phaseOffset=rf_phase/180*pi;
            adc.phaseOffset=rf_phase/180*pi;
            rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
            rf_phase=mod(rf_phase+rf_inc, 360.0);
            %
            seq.addBlock(rf,gz);
            seq.addBlock(gxPre,mr.scaleGrad(gyPre,peScales(i)),gzReph);
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gx,adc);
            %             gyPre.amplitude=-gyPre.amplitude;
            seq.addBlock(mr.makeDelay(delayTR(c)),gxSpoil,mr.scaleGrad(gyPre,-peScales(i)),gzSpoil)
        end
    end
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

%% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness*Nslices]);
seq.setDefinition('Name', 'gre');

seq.write('gre_B0maps_TE1.seq')       % Write to pulseq file

%seq.install('siemens');

%% plot sequence and k-space diagrams

seq.plot('timeRange', [0 5]*TR);

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points

gradSpectrum_pulseq_fct_MPITUB;

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits

rep = seq.testReport;
fprintf([rep{:}]);

