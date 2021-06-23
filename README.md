# MCEGENpiN_radcorr V7b
This V7b version of the newest MC event generator for exclusive single pion electroproduction allows you to generate a massive statistics for the large invariants' scales. Relatively fast and effective procedures in conjunction with model representations make this program a convenient and reliable choice for data analysis in particle physics.

## What's it all about?
One of the biggest parts of any experiments in physics is data analysis. Hadron physics with [CLAS12 spectrometer](https://www.jlab.org/physics/hall-b/clas12) gets pretty tricky when one should deal with its efficiency. This is where programs like this generator come quite handy. Not only do they allow you to restore the original cross-section, but they also can be used as an instrument for event selection development. 

The elaboration of this generator was carried out based on [MAID](https://maid.kph.uni-mainz.de/maid2007/mult.html) representations. As a starting point, we use multipole amplitudes for the charged channels. This data is needed to evaluate the differential cross-section that we further use as weights for event generation.

### Formalism 
Helicity amplitudes were considered the most convenient intermediate stage of the whole data handling. Using the multipole decomposition with Legendre polynomials one can obtain Helicity amplitudes as follows:

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{1}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\theta}\cos{\dfrac{\theta}{2}}\sum_l&space;(E_{l&plus;}&space;-&space;M_{l&plus;}&space;-&space;E_{(l&plus;1)-}&space;-&space;M_{(l&plus;1)-})[P_{l}''(\cos{\theta})&space;-&space;P_{l&plus;1}''(\cos{\theta})]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{1}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\theta}\cos{\dfrac{\theta}{2}}\sum_l&space;(E_{l&plus;}&space;-&space;M_{l&plus;}&space;-&space;E_{(l&plus;1)-}&space;-&space;M_{(l&plus;1)-})[P_{l}''(\cos{\theta})&space;-&space;P_{l&plus;1}''(\cos{\theta})]" title="H_{1} = \dfrac{1}{\sqrt{2}}\sin{\theta}\cos{\dfrac{\theta}{2}}\sum_l (E_{l+} - M_{l+} - E_{(l+1)-} - M_{(l+1)-})[P_{l}''(\cos{\theta}) - P_{l+1}''(\cos{\theta})]" /></a>


<a href="https://www.codecogs.com/eqnedit.php?latex=H_{2}&space;=&space;\dfrac{1}{\sqrt{2}}\cos{\dfrac{\theta}{2}}\sum_l&space;((l&plus;2)E_{l&plus;}&space;-&space;lM_{l&plus;}&space;-&space;(l&plus;2)M_{(l&plus;1)-}&space;&plus;&space;lE_{(l&plus;1)-})[P_{l}'&space;-&space;P_{l&plus;1}']" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{2}&space;=&space;\dfrac{1}{\sqrt{2}}\cos{\dfrac{\theta}{2}}\sum_l&space;((l&plus;2)E_{l&plus;}&space;-&space;lM_{l&plus;}&space;-&space;(l&plus;2)M_{(l&plus;1)-}&space;&plus;&space;lE_{(l&plus;1)-})[P_{l}'&space;-&space;P_{l&plus;1}']" title="H_{2} = \dfrac{1}{\sqrt{2}}\cos{\dfrac{\theta}{2}}\sum_l ((l+2)E_{l+} - lM_{l+} - (l+2)M_{(l+1)-} + lE_{(l+1)-})[P_{l}' - P_{l+1}']" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{3}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\theta}\sin{\dfrac{\theta}{2}}\sum_l&space;(E_{l&plus;}&space;-&space;M_{l&plus;}&space;&plus;&space;E_{(l&plus;1)-}&space;&plus;&space;M_{(l&plus;1)-})[P_{l}''(\cos{\theta})&space;&plus;&space;P_{l&plus;1}''(\cos{\theta})]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{3}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\theta}\sin{\dfrac{\theta}{2}}\sum_l&space;(E_{l&plus;}&space;-&space;M_{l&plus;}&space;&plus;&space;E_{(l&plus;1)-}&space;&plus;&space;M_{(l&plus;1)-})[P_{l}''(\cos{\theta})&space;&plus;&space;P_{l&plus;1}''(\cos{\theta})]" title="H_{3} = \dfrac{1}{\sqrt{2}}\sin{\theta}\sin{\dfrac{\theta}{2}}\sum_l (E_{l+} - M_{l+} + E_{(l+1)-} + M_{(l+1)-})[P_{l}''(\cos{\theta}) + P_{l+1}''(\cos{\theta})]" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{4}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\dfrac{\theta}{2}}\sum_l&space;((l&plus;2)E_{l&plus;}&space;-&space;lM_{l&plus;}&space;&plus;&space;(l&plus;2)M_{(l&plus;1)-}&space;-&space;lE_{(l&plus;1)-})[P_{l}'&space;&plus;&space;P_{l&plus;1}']" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{4}&space;=&space;\dfrac{1}{\sqrt{2}}\sin{\dfrac{\theta}{2}}\sum_l&space;((l&plus;2)E_{l&plus;}&space;-&space;lM_{l&plus;}&space;&plus;&space;(l&plus;2)M_{(l&plus;1)-}&space;-&space;lE_{(l&plus;1)-})[P_{l}'&space;&plus;&space;P_{l&plus;1}']" title="H_{4} = \dfrac{1}{\sqrt{2}}\sin{\dfrac{\theta}{2}}\sum_l ((l+2)E_{l+} - lM_{l+} + (l+2)M_{(l+1)-} - lE_{(l+1)-})[P_{l}' + P_{l+1}']" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{5}&space;=\dfrac{Q}{|k|}&space;\cos(\dfrac{\theta}{2})\sum_l&space;(l&plus;1)(S_{l&plus;}&space;&plus;&space;S_{(l&plus;1)-})[P_{l}'(\cos(\theta))&space;-&space;P_{l&plus;1}'(\cos(\theta))]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{5}&space;=\dfrac{Q}{|k|}&space;\cos(\dfrac{\theta}{2})\sum_l&space;(l&plus;1)(S_{l&plus;}&space;&plus;&space;S_{(l&plus;1)-})[P_{l}'(\cos(\theta))&space;-&space;P_{l&plus;1}'(\cos(\theta))]" title="H_{5} =\dfrac{Q}{|k|} \cos(\dfrac{\theta}{2})\sum_l (l+1)(S_{l+} + S_{(l+1)-})[P_{l}'(\cos(\theta)) - P_{l+1}'(\cos(\theta))]" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=H_{6}&space;=\dfrac{Q}{|k|}&space;\sin(\dfrac{\theta}{2})\sum_l&space;(l&plus;1)(S_{l&plus;}&space;-&space;S_{(l&plus;1)-})[P_{l}'(\cos(\theta))&space;&plus;&space;P_{l&plus;1}'(\cos(\theta))]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?H_{6}&space;=\dfrac{Q}{|k|}&space;\sin(\dfrac{\theta}{2})\sum_l&space;(l&plus;1)(S_{l&plus;}&space;-&space;S_{(l&plus;1)-})[P_{l}'(\cos(\theta))&space;&plus;&space;P_{l&plus;1}'(\cos(\theta))]" title="H_{6} =\dfrac{Q}{|k|} \sin(\dfrac{\theta}{2})\sum_l (l+1)(S_{l+} - S_{(l+1)-})[P_{l}'(\cos(\theta)) + P_{l+1}'(\cos(\theta))]" /></a>

These amplitudes are further used for the structure functions evaluation, which one requires for the differential cross-section calculation.

<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_{T}&space;=&space;\dfrac{q}{2K}(|H_1|^2&space;&plus;&space;|H_2|^2&space;&plus;&space;|H_3|^2&space;&plus;&space;|H_4|^2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_{T}&space;=&space;\dfrac{q}{2K}(|H_1|^2&space;&plus;&space;|H_2|^2&space;&plus;&space;|H_3|^2&space;&plus;&space;|H_4|^2)" title="\sigma_{T} = \dfrac{q}{2K}(|H_1|^2 + |H_2|^2 + |H_3|^2 + |H_4|^2)" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_{L}&space;=&space;\dfrac{q}{K}(|H_5|^2&space;&plus;&space;|H_6|^2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_{L}&space;=&space;\dfrac{q}{K}(|H_5|^2&space;&plus;&space;|H_6|^2)" title="\sigma_{L} = \dfrac{q}{K}(|H_5|^2 + |H_6|^2)" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_{TT}&space;=&space;\dfrac{q}{K}Re(H_3H_2^*&space;-&space;H_4H_1^*)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_{TT}&space;=&space;\dfrac{q}{K}Re(H_3H_2^*&space;-&space;H_4H_1^*)" title="\sigma_{TT} = \dfrac{q}{K}Re(H_3H_2^* - H_4H_1^*)" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_{LT}&space;=&space;-\dfrac{q}{\sqrt{2}K}Re((H_1&space;-&space;H_4)H_5^*&space;&plus;&space;(H_2&space;&plus;&space;H_3)H_6^*)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_{LT}&space;=&space;-\dfrac{q}{\sqrt{2}K}Re((H_1&space;-&space;H_4)H_5^*&space;&plus;&space;(H_2&space;&plus;&space;H_3)H_6^*)" title="\sigma_{LT} = -\dfrac{q}{\sqrt{2}K}Re((H_1 - H_4)H_5^* + (H_2 + H_3)H_6^*)" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_{LT'}&space;=&space;-\dfrac{q}{\sqrt{2}K}Im((H_1&space;-&space;H_4)H_5^*&space;&plus;&space;(H_2&space;&plus;&space;H_3)H_6^*)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_{LT'}&space;=&space;-\dfrac{q}{\sqrt{2}K}Im((H_1&space;-&space;H_4)H_5^*&space;&plus;&space;(H_2&space;&plus;&space;H_3)H_6^*)" title="\sigma_{LT'} = -\dfrac{q}{\sqrt{2}K}Im((H_1 - H_4)H_5^* + (H_2 + H_3)H_6^*)" /></a>

where K and k and q are respectively, the photon equivalent energy and the virtual photon and pion 3-momenta in γN c.m.s. For unpolarized particles and for a longitudinally polarized electron beam, the φ-dependence of the γN → Nπ cross section can be specified in the following way:

<a href="https://www.codecogs.com/eqnedit.php?latex=\dfrac{d\sigma}{d\Omega}&space;=&space;\sigma_{T}&space;&plus;&space;\epsilon&space;\sigma_{L}&space;&plus;&space;\epsilon&space;\sigma_{TT}\cos{2\varphi}&space;&plus;&space;\sqrt{2\epsilon&space;(1&space;&plus;&space;\epsilon)}\sigma_{LT}\cos{\varphi}&space;&plus;&space;h\sqrt{2\epsilon(1&space;-&space;\epsilon)}\sigma_{LT'}\sin{\varphi}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dfrac{d\sigma}{d\Omega}&space;=&space;\sigma_{T}&space;&plus;&space;\epsilon&space;\sigma_{L}&space;&plus;&space;\epsilon&space;\sigma_{TT}\cos{2\varphi}&space;&plus;&space;\sqrt{2\epsilon&space;(1&space;&plus;&space;\epsilon)}\sigma_{LT}\cos{\varphi}&space;&plus;&space;h\sqrt{2\epsilon(1&space;-&space;\epsilon)}\sigma_{LT'}\sin{\varphi}" title="\dfrac{d\sigma}{d\Omega} = \sigma_{T} + \epsilon \sigma_{L} + \epsilon \sigma_{TT}\cos{2\varphi} + \sqrt{2\epsilon (1 + \epsilon)}\sigma_{LT}\cos{\varphi} + h\sqrt{2\epsilon(1 - \epsilon)}\sigma_{LT'}\sin{\varphi}" /></a>

The differential cross section of the electroproduction of pions off nucleons in the one-photon approximation was used as a weight for each specific event.

<a href="https://www.codecogs.com/eqnedit.php?latex=\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega}&space;=&space;\Gamma&space;\dfrac{d\sigma}{d\Omega}&space;\;\;\;\;\;\;\;\;&space;\Gamma&space;=&space;\dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2&space;-&space;m^2)E_f}{2mE_i}\dfrac{1}{1&space;-&space;\epsilon}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega}&space;=&space;\Gamma&space;\dfrac{d\sigma}{d\Omega}&space;\;\;\;\;\;\;\;\;&space;\Gamma&space;=&space;\dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2&space;-&space;m^2)E_f}{2mE_i}\dfrac{1}{1&space;-&space;\epsilon}" title="\dfrac{d\sigma}{dE_fd\Omega_{e}d\Omega} = \Gamma \dfrac{d\sigma}{d\Omega} \;\;\;\;\;\;\;\; \Gamma = \dfrac{\alpha}{2\pi^2Q^2}\dfrac{(W^2 - m^2)E_f}{2mE_i}\dfrac{1}{1 - \epsilon}" /></a>

### Radiative Corrections

For RE(radiative effects) simulations the [Mo and Tsai](https://inspirehep.net/literature/52657) approach was chosen. Using the peak approximation we were able to implement the following corrections:
  1. Weights re-evaluation according to RC
  2. The radiative tail simulation
  3. Outgoing rad. photon generation

### Main idea

The general procedure for the event generation contains 2 parts:
* Weight calculation
* Particles kinematics

In the first case, the program decides whether it's the soft or hard region for the RC and completes the differential cross section calculation.
The additional factor for the weight is also calculated here. The program completes the event generation cycle with the kinematics evaluation. Thus it writes all the information about the specific event in the output file using the ["Lund"](https://gemc.jlab.org/gemc/html/documentation/generator/lund.html) format. 

As a result, for each <a href="https://www.codecogs.com/eqnedit.php?latex=\{W,&space;Q^2,&space;cos(\theta),&space;\varphi\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\{W,&space;Q^2,&space;cos(\theta),&space;\varphi\}" title="\{W, Q^2, cos(\theta), \varphi\}" /></a> point and <a href="https://www.codecogs.com/eqnedit.php?latex=E_{beam},&space;h" target="_blank"><img src="https://latex.codecogs.com/gif.latex?E_{beam},&space;h" title="E_{beam}, h" /></a> (beam polarization) values, this program creates a bunch of outgoing particles with the appropriate kinematic values.

## Usage

## Some histograms

<p align="left"> <img src="https://komarev.com/ghpvc/?username=maksaska&label=Profile%20views&color=0e75b6&style=flat" alt="maksaska" /> <img src="https://img.shields.io/badge/MSU-SINP-blue" /> <img src="https://img.shields.io/badge/JLab-red" /> </p>
