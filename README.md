# MCEGENpiN_radcorr V7b
This V7b version of the newest MC event generator for exclusive single pion electroproduction allows you to generate a massive statistics for the large invariants' scales. Relatively fast and effective procedures in conjunction with model representations make this program a convenient and reliable choice for data analysis in particle physics.

## What's it all about?
One of the biggest parts of any experiments in physics is data analysis. Hadron physics with CLAS12 spectrometer gets pretty tricky when one should deal with its efficiency. This is where programs like this generator come quite handy. Not only do they allow you to restore the original cross-section, but they also can be used as an instrument for event selection development. 

The elaboration of this generator was carried out based on MAID representations. As a starting point, we use multipole amplitudes for the charged channels. This data is needed to evaluate the differential cross-section that we further use as weights for event generation.

### Formalism 
Helicity amplitudes were considered the most convenient intermediate stage of the whole data handling. Using the multipole decomposition with Legendre polynomials one can obtain Helicity amplitudes as follows:

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{1}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\theta}\cos{\dfrac{\theta}{2}}\sum_l&space;(E_{l&plus;}&space;-&space;M_{l&plus;}&space;-&space;E_{(l&plus;1)-}&space;-&space;M_{(l&plus;1)-})[P_{l}''(\cos{\theta})&space;-&space;P_{l&plus;1}''(\cos{\theta})]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{1}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\theta}\cos{\dfrac{\theta}{2}}\sum_l&space;(E_{l&plus;}&space;-&space;M_{l&plus;}&space;-&space;E_{(l&plus;1)-}&space;-&space;M_{(l&plus;1)-})[P_{l}''(\cos{\theta})&space;-&space;P_{l&plus;1}''(\cos{\theta})]" title="H_{1} = \dfrac{1}{\sqrt{2}}\sin{\theta}\cos{\dfrac{\theta}{2}}\sum_l (E_{l+} - M_{l+} - E_{(l+1)-} - M_{(l+1)-})[P_{l}''(\cos{\theta}) - P_{l+1}''(\cos{\theta})]" /></a>


<a href="https://www.codecogs.com/eqnedit.php?latex=H_{2}&space;=&space;\dfrac{1}{\sqrt{2}}\cos{\dfrac{\theta}{2}}\sum_l&space;((l&plus;2)E_{l&plus;}&space;-&space;lM_{l&plus;}&space;-&space;(l&plus;2)M_{(l&plus;1)-}&space;&plus;&space;lE_{(l&plus;1)-})[P_{l}'&space;-&space;P_{l&plus;1}']" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{2}&space;=&space;\dfrac{1}{\sqrt{2}}\cos{\dfrac{\theta}{2}}\sum_l&space;((l&plus;2)E_{l&plus;}&space;-&space;lM_{l&plus;}&space;-&space;(l&plus;2)M_{(l&plus;1)-}&space;&plus;&space;lE_{(l&plus;1)-})[P_{l}'&space;-&space;P_{l&plus;1}']" title="H_{2} = \dfrac{1}{\sqrt{2}}\cos{\dfrac{\theta}{2}}\sum_l ((l+2)E_{l+} - lM_{l+} - (l+2)M_{(l+1)-} + lE_{(l+1)-})[P_{l}' - P_{l+1}']" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{3}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\theta}\sin{\dfrac{\theta}{2}}\sum_l&space;(E_{l&plus;}&space;-&space;M_{l&plus;}&space;&plus;&space;E_{(l&plus;1)-}&space;&plus;&space;M_{(l&plus;1)-})[P_{l}''(\cos{\theta})&space;&plus;&space;P_{l&plus;1}''(\cos{\theta})]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{3}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\theta}\sin{\dfrac{\theta}{2}}\sum_l&space;(E_{l&plus;}&space;-&space;M_{l&plus;}&space;&plus;&space;E_{(l&plus;1)-}&space;&plus;&space;M_{(l&plus;1)-})[P_{l}''(\cos{\theta})&space;&plus;&space;P_{l&plus;1}''(\cos{\theta})]" title="H_{3} = \dfrac{1}{\sqrt{2}}\sin{\theta}\sin{\dfrac{\theta}{2}}\sum_l (E_{l+} - M_{l+} + E_{(l+1)-} + M_{(l+1)-})[P_{l}''(\cos{\theta}) + P_{l+1}''(\cos{\theta})]" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{4}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\dfrac{\theta}{2}}\sum_l&space;((l&plus;2)E_{l&plus;}&space;-&space;lM_{l&plus;}&space;&plus;&space;(l&plus;2)M_{(l&plus;1)-}&space;-&space;lE_{(l&plus;1)-})[P_{l}'&space;&plus;&space;P_{l&plus;1}']" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{4}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\dfrac{\theta}{2}}\sum_l&space;((l&plus;2)E_{l&plus;}&space;-&space;lM_{l&plus;}&space;&plus;&space;(l&plus;2)M_{(l&plus;1)-}&space;-&space;lE_{(l&plus;1)-})[P_{l}'&space;&plus;&space;P_{l&plus;1}']" title="H_{4} = \dfrac{1}{\sqrt{2}}\sin{\dfrac{\theta}{2}}\sum_l ((l+2)E_{l+} - lM_{l+} + (l+2)M_{(l+1)-} - lE_{(l+1)-})[P_{l}' + P_{l+1}']" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{5}&space;=\dfrac{Q}{|k|}&space;\cos(\dfrac{\theta}{2})\sum_l&space;(l&plus;1)(S_{l&plus;}&space;&plus;&space;S_{(l&plus;1)-})[P_{l}'(\cos(\theta))&space;-&space;P_{l&plus;1}'(\cos(\theta))]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{5}&space;=\dfrac{Q}{|k|}&space;\cos(\dfrac{\theta}{2})\sum_l&space;(l&plus;1)(S_{l&plus;}&space;&plus;&space;S_{(l&plus;1)-})[P_{l}'(\cos(\theta))&space;-&space;P_{l&plus;1}'(\cos(\theta))]" title="H_{5} =\dfrac{Q}{|k|} \cos(\dfrac{\theta}{2})\sum_l (l+1)(S_{l+} + S_{(l+1)-})[P_{l}'(\cos(\theta)) - P_{l+1}'(\cos(\theta))]" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{6}&space;=\dfrac{Q}{|k|}&space;\sin(\dfrac{\theta}{2})\sum_l&space;(l&plus;1)(S_{l&plus;}&space;-&space;S_{(l&plus;1)-})[P_{l}'(\cos(\theta))&space;&plus;&space;P_{l&plus;1}'(\cos(\theta))]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{6}&space;=\dfrac{Q}{|k|}&space;\sin(\dfrac{\theta}{2})\sum_l&space;(l&plus;1)(S_{l&plus;}&space;-&space;S_{(l&plus;1)-})[P_{l}'(\cos(\theta))&space;&plus;&space;P_{l&plus;1}'(\cos(\theta))]" title="H_{6} =\dfrac{Q}{|k|} \sin(\dfrac{\theta}{2})\sum_l (l+1)(S_{l+} - S_{(l+1)-})[P_{l}'(\cos(\theta)) + P_{l+1}'(\cos(\theta))]" /></a>

### Main idea

## Usage

## Some histograms

<p align="left"> <img src="https://komarev.com/ghpvc/?username=maksaska&label=Profile%20views&color=0e75b6&style=flat" alt="maksaska" /> <img src="https://img.shields.io/badge/MSU-SINP-blue" /> <img src="https://img.shields.io/badge/JLab-red" /> </p>
