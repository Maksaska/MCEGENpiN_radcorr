# MCEGENpiN_radcorr V7b
This V7b version of the newest MC event generator for exclusive single pion electroproduction allows you to generate a massive statistics for the large invariants' scales. Relatively fast and effective procedures in conjunction with model representations make this program a convenient and reliable choice for data analysis in particle physics.

## What's it all about?
One of the biggest parts of any experiments in physics is data analysis. Hadron physics with CLAS12 spectrometer gets pretty tricky when one should deal with its efficiency. This is where programs like this generator come quite handy. Not only do they allow you to restore the original cross-section, but they also can be used as an instrument for event selection development. 

The elaboration of this generator was carried out based on MAID representations. As a starting point, we use multipole amplitudes for the charged channels. This data is needed to evaluate the differential cross-section that we further use as weights for event generation.

### Formalism 
Helicity amplitudes were considered the most convenient intermediate stage of the whole data handling. Using the multipole decomposition with Legendre polynomials one can obtain Helicity amplitudes as follows:

<img src="https://bit.ly/3gAuSpx" align="center" border="0" alt="H_{1} = \dfrac{1}{\sqrt{2}}\sin{\theta}\cos{\dfrac{\theta}{2}}\sum_l (E_{l+} - M_{l+} - E_{(l+1)-} - M_{(l+1)-})[P_{l}''(\cos{\theta}) - P_{l+1}''(\cos{\theta})]" width="635" height="50" />
<img src="https://bit.ly/3q7lhcP" align="center" border="0" alt="H_{2} = \dfrac{1}{\sqrt{2}}\cos{\dfrac{\theta}{2}}\sum_l ((l+2)E_{l+} - lM_{l+} - (l+2)M_{(l+1)-} + lE_{(l+1)-})[P_{l}'(\cos{\theta}) - P_{l+1}'(\cos{\theta})]" width="700" height="50" /
<img src="https://bit.ly/3wBuSed" align="center" border="0" alt="H_{3} = \dfrac{1}{\sqrt{2}}\sin{\theta}\sin{\dfrac{\theta}{2}}\sum_l (E_{l+} - M_{l+} + E_{(l+1)-} + M_{(l+1)-})[P_{l}''(\cos{\theta}) + P_{l+1}''(\cos{\theta})]" width="631" height="50" />
<img src="https://bit.ly/3iVKKEF" align="center" border="0" alt="H_{4} = \dfrac{1}{\sqrt{2}}\sin{\dfrac{\theta}{2}}\sum_l ((l+2)E_{l+} - lM_{l+} + (l+2)M_{(l+1)-} - lE_{(l+1)-})[P_{l}'(\cos{\theta}) + P_{l+1}'(\cos{\theta})]" width="696" height="50" />
<img src="https://bit.ly/3iPwugl" align="center" border="0" alt="H_{5} =\dfrac{Q}{|k|}  \cos(\dfrac{\theta}{2})\sum_l (l+1)(S_{l+} + S_{(l+1)-})[P_{l}'(\cos(\theta)) - P_{l+1}'(\cos(\theta))]" width="543" height="50" />
<img src="https://bit.ly/3zCkPHK" align="center" border="0" alt="H_{6} =\dfrac{Q}{|k|}  \sin(\dfrac{\theta}{2})\sum_l (l+1)(S_{l+} - S_{(l+1)-})[P_{l}'(\cos(\theta)) + P_{l+1}'(\cos(\theta))]" width="539" height="50" />

### Main idea

## Usage

## Some histograms


<p align="left"> <img src="https://komarev.com/ghpvc/?username=maksaska&label=Profile%20views&color=0e75b6&style=flat" alt="maksaska" /> <img src="https://img.shields.io/badge/MSU-SINP-blue" /> <img src="https://img.shields.io/badge/JLab-red" /> </p>
