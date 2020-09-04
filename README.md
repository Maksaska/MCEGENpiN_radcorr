The purpose of this program is to generate the events for Pi0p an Pi+n channels of meson electroproduction reaction in the kinematic range of W:[1.08, 2] GeV and Q^2: [0, 10] GeV^2. For Q^2: [5, 10] GeV^2 the data was extrapolated. For the extrapolation the following procedure was used : All multipoles for W: 1.08 - 2 GeV and Q2: 5 - 10 GeV^2 were copied from Q^2  = 5 GeV^2 and corresponding W with sqrt(5/Q2)^(n/2) factor.
For input data you can use the "input.txt" file. 

To launch the generator use "Run" or "Start" file. "Run" = "Start" + compilation 

For input file: 

Enter the information you are interested in in the following order: Beam energy (GeV), left border of the range W (GeV), right border of the range W (GeV), left border of the range Q^2 (GeV^2), right border of the range Q^2 (GeV), beam polarization h(+-1, 0), n factor for the extrapolation, Number of genereted events, Histogramm request(1 - Yes; 0 - No), Total amount of ".lund" files, Channel (2 - Pi+n; 1 - Pi0p + Pi0 decay; 0 - pi0p), Interpolation mode (0 - linear; 1 - quad), "weights" mode (1 - Flat distribution over the entire kinematic range. Weights(cross-sections = [mubr]) are being recorded in .lund file; 0 - Events are recorded taking into account the peculiarities of their distribution in the kinetic region). Enter all input data separated by a space.

Example: 10.6 1.08 2 0 10 0 4 100000 1 1 2 0 1

This input stands for 10.6 GeV beam energy(beam polarization is h = 0) with the kinematic range of W:[1.08, 2] GeV and Q^2: [0, 10] GeV^2 for 100000 events, that will be recorded in 1 file "pin_1.08_2_0_5_(1).lund". Additional 2d histogram (W, Q^2) will be created. Cross-section for random bin (W, Q^2) will be calculated with the help of linear interpolation. If the extrapolation is necessary the n = 4 factor will be used. Flat distribution over the entire kinematic range. Weights are recorded in .lund file.

