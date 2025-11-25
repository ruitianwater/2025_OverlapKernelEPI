%% Fucntion to plot gradient spectrum along with known acoustic resonances
% acoustic resonances are obatained from imprint.exe of IDEA
% gradient spectrum calculation is from pulseq matlab repo
% last modified: 2024-07-03(praveen)
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.

%% PULSEQ sequences

%check pulseq on path
assert(exist('gradSpectrum.m','file'),'add pulseq(v1.4) (https://github.com/pulseq/pulseq.git) to path');

system = mr.opts('rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
                 'adcDeadTime', 20e-6);

seq=mr.Sequence(system);              % Create a new sequence object
[fn,pn] = uigetfile('*.seq','please select your *.seq file','MultiSelect','off');
seq.read(fullfile([pn,fn]));

gradSpectrum_pulseq_fct(seq,fn);

return;
%% SIEMENS sequences
% 
%check dsv_read on path
assert(exist('dsv_read.m','file'),'add dsv_read.m to path');

[pn_dsv] = uigetdir('please selct directory that contains the *.dsv files');
filename{1} = [pn_dsv,'\SimulationProtocol_GRX.dsv'];
filename{2} = [pn_dsv,'\SimulationProtocol_GRY.dsv'];
filename{3} = [pn_dsv,'\SimulationProtocol_GRZ.dsv'];
fprintf('reading data from dsv ...');
for nn=1:3
    [data{nn}(2,:), header{nn}]= dsv_read(filename{nn});
    data{nn}(1,:) = linspace(0,header{nn}.HORIDELTA*10^-6 * size(data{nn},2), size(data{nn},2)   );
end
fprintf(' ok done!\n');
%some dummy variables to be able to call the function 
seq.sys.gradRasterTime = header{nn}.HORIDELTA * 10^-6;
mod_pn = strfind(pn_dsv,'\');
mod_pn = pn_dsv(mod_pn(end)+1:end);
fn = sprintf('FromSiemens_%s.dsv',mod_pn);

fprintf('starting gradient calculations ...\n');
gradSpectrum_pulseq_fct(seq,fn,data);


%% % NESTED FUNCTIONS ONLY %%%

function gradSpectrum_pulseq_fct(seq,fn,varargin)
% this (ab)uses the gradSpectrum.m batch of pulseq 

% example script to plot gradient frequency spectrum
% expects: 
%  "seq" object to be already populated 
%  "sys" object to contain system specs
% furthermore, on Siemens you need the *.asc file for your gradient system

try
    dt=sys.gradRasterTime; % time raster
catch
    dt=seq.sys.gradRasterTime;
    fprintf('WARNING:\n\treading sys properties from your seq file!\n\traster time is assumed to be %.8f [sec?]\n',dt);
end
fmax=10000/2; %10kHz
nwin=5000; % 0.05s
os=3; % frequency oversampling for prettier peaks
ascName=[]; % this disables the display of the system's resonance frequences
% ascName='idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'; % 3T prisma
% ascName='idea/asc/MP_GradSys_P034_X60.asc'; % 3T cima.X

if ischar(ascName)
    ascData=mr.Siemens.readasc(ascName);
end

faxis=(0:(nwin/2-1))/nwin/dt/os;
nfmax=sum(faxis<=fmax);

if nargin < 3
    wave_data=seq.waveforms_and_times();
else
    fprintf('*** WARNING ***\n\tNOT reading from a *.seq but from an external file!\n');
    wave_data = varargin{1};
end

tmax=max([wave_data{1}(1,end) wave_data{2}(1,end) wave_data{3}(1,end)]);
nt=ceil(tmax/dt);
tmax=nt*dt;

gw=zeros(3,nt);
for i=1:3
    gw(i,:)=interp1(wave_data{i}(1,:),wave_data{i}(2,:),((1:nt)-0.5)*dt,'linear',0);
    % alternative (to be checked in the future)
    % it is actually much more appropriate to calculate the spectrium of
    % the derivative(!) of the gradient wave form and not the waveform
    % itself, at least for the cound/noise of the gradients...
    %gw(i,1:end-1)=diff(interp1(wave_data{i}(1,:),wave_data{i}(2,:),((1:nt)-0.5)*dt,'linear',0));
end

gs=[];

ng=size(gw,1);
for g=1:ng
    x=gw(g,:);
    nx = length(x);

    nx=ceil(nx/nwin)*nwin;
    if nx>length(x) 
        x=[x, zeros(1,nx-length(x))]; % zerofill
    end

    nseg1=nx/nwin; 
    xseg=zeros(nseg1*2-1,nwin*os); 

    xseg(1:2:end,1:nwin)=reshape(x,[nwin,nseg1])';
    if nseg1>1
        xseg(2:2:end,1:nwin)=reshape(x(1+nwin/2:end-nwin/2),[nwin,nseg1-1])';
    end

    xseg_dc=mean(xseg,2);
    xseg=xseg-xseg_dc(:,ones(1,nwin*os));

    if nseg1>1 % WARNING: thisintroduces inconsistency between short and long sequences in term os the peak amplitudes
        cwin=0.5*(1-cos(2*pi*(1:nwin)/nwin));
        xseg(:,1:nwin)=xseg(:,1:nwin).*cwin(ones(size(xseg,1),1),:);
    end

    fseg=abs(fft(xseg,[],2));
    fseg=fseg(:,1:end/2); 
        
    if nseg1>1 
        gs = [gs; mean(fseg.^2).^0.5]; % sos
        %figure; plot(faxis(1:nfmax),sum(fseg(:,1:nfmax).^2).^0.5);
    else
        gs = [gs; abs(fseg)]; % add abs
    end    
end

figure('Units','normalized','OuterPosition',[0 0 1 1]); 
plot(faxis(1:nfmax),gs(:,1:nfmax),'LineWidth',3);
%%figure; plot(faxis(1:nfmax),sum(gs(:,1:nfmax))); % abs-sum
hold on; plot(faxis(1:nfmax),sum(gs(:,1:nfmax).^2).^0.5,'LineWidth',3); % sos
xlabel('frequency / Hz');

if ischar(ascName)
    hold on; 
    for i=1:length(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency)
        if ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)>0
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i),'-');
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)-ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i)/2,'--');
            xline(ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceFrequency(i)+ascData.asGPAParameters(1).sGCParameters.aflAcousticResonanceBandwidth(i)/2,'--');
        end
    end
end

legend({'Gx','Gy','Gz','Gtot'});


%nov = floor(nsc/2);
%nff = max(256,2^nextpow2(nsc));
%t = spectrogram(x,hamming(nsc),nov);%,nff);
%t = spectrogram(x,rectwin(nsc),nov);
%maxerr = max(abs(abs(t(:))-abs(s(:))))

%smueller: add system limits as overlay, code initially written by pvalal
System2Plot = questdlg('which system you want?','GUI','Prisma','Prisma');

% acoustic resonances
switch System2Plot
    case 'Prisma'
        f=[590 1140];% Hz
        BW=[100 220]; %Hz
end
x=[f-0.5*BW;f-0.5*BW;f+0.5*BW;f+0.5*BW];
tax=gca;
y = [zeros(1,numel(f)); ones(1,numel(f))*tax.YLim(2); ones(1,numel(f))*tax.YLim(2); zeros(1,numel(f));];
c = [ones(1,numel(f))*-1;ones(1,numel(f))*-1;ones(1,numel(f)); ones(1,numel(f))];
hold on,
patch(x,y,c,'FaceAlpha',0.25)
colormap(vertcat(jet,flipud(jet)))
tax.Title.String = sprintf('System is %s for sequence "%s"',System2Plot,strrep(fn,'_','\_'));
tax.FontSize=20;
grid on;grid minor;
set(gcf,'color','white')
legend({'Gx','Gy','Gz','Gtot','DANGER'});
%%%%%
end % main function