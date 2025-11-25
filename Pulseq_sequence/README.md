==========================================================
Sequence tests

Run the following MATLAB files for generating sequences


rtian\_00\_ACS\_prescan.m (GRE, ACS reference scan)
rtian\_00\_B0maps\_TE1.m (GRE, B0 mapping TE1)
rtian\_00\_B0maps\_TE2.m (GRE, B0 mapping TE2)
rtian\_01\_ss\_SE\_epi\_b1k.m (single-shot EPI)
rtian\_02\_2ps\_SE\_epi\_b1k.m (2-shot mosaic EPI)
rtian\_03\_4mo\_SE\_epi\_b1k.m (4-shot mosaic EPI)
rtian\_04\_5rs\_SE\_epi\_b1k.m (5-shot readout-segmented EPI)
rtian\_05\_5pi\_SE\_epi\_b1k.m (5-shot phase-interleaved EPI)



Rename the files (xxx.seq), run these sequences in the MR scanners with Pulseq installed.



To reconstruct the data acquired by these sequences, please save all the data (Twix in Siemens) and .seq files properly, and check the scripts in the "Pulseq\_recon" folder.



The pulseq-master folder contains the Pulseq package that was used to generate the multi-shot EPI sequences and to facilitate image reconstruction. You could download the lastest version of Pulseq online for replacement.

========================================================

Please cite the following conference/workshop abstracts, and the incoming research articles for this work:



Tian R, Uecker M, Zaitsev M, Scheffler K. Overlap-kernel EPI: estimating MRI shot-to-shot phase variations by shifted-kernel extraction from overlapped regions at arbitrary k-space locations. Magn Reson Med, 2025. doi: 10.1002/mrm.70196


Tian R, Uecker M, Scheffler K. Navigator-free multi-shot EPI with shift-invariant kernel extraction in subspace. Power pitch. Program #1081. In proceedings of the annual meeting of ISMRM, Singapore, 2024.


Tian R, Uecker M, Zaitsev M, Scheffler K. Pulseq Implementation of Overlapped Readout Segmented EPI with Phase Fluctuations Corrected by Shift-Invariant Kernel Extraction. Traditional poster #18. ISMRM Workshop on 40 Years of Diffusion: Past, Present \& Future Perspectives, Kyoto, Japan, 2025.

========================================================
Modified from the experiment 2025-11-06

