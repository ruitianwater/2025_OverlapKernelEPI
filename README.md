This is a package of example codes and data for the Overlap-Kernel EPI paper, including Pulseq sequence and offline reconstruction, under BSD 3-Clause License.

The Overlap-Kernel EPI is a multi-shot EPI technique with inter-shot k-space overlaps for magnetic resonance imaging. 
These overlap regions serve as calibration region, where a GRAPPA/ESPIRiT-type operation can be used to explicitly extract the notorious shot-to-shot phase fluctuation kernels/maps, which are then incorporated into multi-shot reconstruction to avoid artifacts. 
This allows self-navigation for various multi-shot trajectories, without separate 2D navigators given sufficient SNR. 
This package contains options for different EPI trajectories, with and without 2D navigators, spin-echo diffusion-weighting or gradient-echo acquisitions, different shot-to-shot phase corrections, and many other features for sequence design and reconstruction.

Please download the .zip file, unzip it, and check the multi-shot EPI sequence generation and the related offline reconstruction. A Pulseq package is also included, which was used to generate sequences and support the reconstruction.

Please cite the following paper:

Tian R, Uecker M, Zaitsev M, Scheffler K. Overlap-kernel EPI: estimating MRI shot-to-shot phase variations by shifted-kernel extraction from overlapped regions at arbitrary k-space locations. Magn Reson Med, 2025. doi: 10.1002/mrm.70196

For further questions, please contact Rui Tian, rui.tian@tuebingen.mpg.de/ruitianwater@outlook.com
