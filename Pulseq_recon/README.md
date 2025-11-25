==========================================================
Reconstruction tests

Run the following MATLAB files for testing reconstruction



single-shot EPI: Script1\_ss\_EPI\_Recon.m



2-shot phase mosaic segmented EPI: Script2\_2ps\_EPI\_Recon.m



4-shot mosaic segmented EPI: Script3\_4mo\_EPI\_Recon.m



5-shot readout segmented EPI: Script4\_5rs\_EPI\_Recon.m



6-shot phase interleaved EPI: Script5\_5pi\_EPI\_Recon.m



For multi-shot EPI, a few options can be changed following the code comments, for testing purposes:


para.bIntSCorr % To use self-navigated inter-shot corrections or not


para.type\_subspace = 1; % The self-navigated correction algorithms


para.bNav = 0; % To use 2D navigator inter-shot corrections or not



The pulseq-master folder contains the Pulseq package that was used to generate the multi-shot EPI sequences and to facilitate image reconstruction. You could download the lastest version of Pulseq online for replacement.



========================================================



Please cite the following conference/workshop abstracts, and the incoming research articles for this work:



Tian R, Uecker M, Zaitsev M, Scheffler K. Overlap-kernel EPI: estimating MRI shot-to-shot phase variations by shifted-kernel extraction from overlapped regions at arbitrary k-space locations. Magn Reson Med, 2025. doi: 10.1002/mrm.70196


Tian R, Uecker M, Scheffler K. Navigator-free multi-shot EPI with shift-invariant kernel extraction in subspace. Power pitch. Program #1081. In proceedings of the annual meeting of ISMRM, Singapore, 2024.


Tian R, Uecker M, Zaitsev M, Scheffler K. Pulseq Implementation of Overlapped Readout Segmented EPI with Phase Fluctuations Corrected by Shift-Invariant Kernel Extraction. Traditional poster #18. ISMRM Workshop on 40 Years of Diffusion: Past, Present \& Future Perspectives, Kyoto, Japan, 2025.

========================================================
Modified from the experiment 2025-11-06



