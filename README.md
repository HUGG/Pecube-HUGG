Pecube-HUGG
===========

This is a repository for the Helsinki University Geodynamics Group (HUGG) version of the thermokinematic and thermochronometer age prediction program Pecube (Pecube-HUGG).
The HUGG version of Pecube was forked from Pecube version 3.0, distributed by Jean Braun in September 2010.

Licensing
---------

Pecube is open source, but attribution is requested by citing the articles below.

Attribution
-----------

When publishing or presenting results from Pecube experiments please cite the following article(s):

* Braun, J. (2003). Pecube: A new finite element code to solve the heat transport equation in three dimensions in the Earth's crust including the effects of a time-varying, finite amplitude surface topography. _Compututers and Geosciences_, 29:787-794.
* Braun, J., der Beek, van, P., Valla, P., Robert, X., Herman, F., Glotzbach, C., et al. (2012). Quantifying rates of landscape evolution and tectonic processes by thermochronology and numerical modeling of crustal heat transport using PECUBE. _Tectonophysics_, 524-525(0), 1-28. doi:10.1016/j.tecto.2011.12.035

Example model and data
----------------------

This version of Pecube includes example model input files ([fault_parameters.txt](input/fault_parameters.txt), [topo_parameters.txt](input/topo_parameters.txt)), and topographic and bedrock age data from western Bhutan.
This data and the model designs are described in detail in Coutand et al. (2014).

* Coutand, I., Whipp, D. M., Grujic, D., Bernet, M., Fellin, M. G., Bookhagen, B., Landry, K., Ghalley, S. & Duncan, C. (2014). [Geometry and kinematics of the Main Himalayan Thrust and Neogene crustal exhumation in the Bhutanese Himalaya derived from inversion of multithermochronologic data](https://dx.doi.org/10.1002/2013JB010891). _Journal of Geophysical Research: Solid Earth_, 119(2), 1446-1481.

Data inversion using the Neighbourhood Algorithm
------------------------------------------------

If performing inversions of thermochronometer data using the [Neighbourhood Algorithm](http://rses.anu.edu.au/~malcolm/na/), be aware that as described in the software package, use of the Neighbourhood Algorithm is restricted to non-commercial research and teaching.
Any commercial use requires explicit written permission.
More information about the Neighbourhood Algorithm and its use can be found on the [Neighbourhood Algorithm web page](http://rses.anu.edu.au/~malcolm/na/).
