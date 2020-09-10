This folder will contain the data, programs and pdf for the paper:<br>
Quénot, Georges & Pakleza, Jaroslaw & Kowalewski, Tomasz. (1998). Particle Image Velocimetry with Optical Flow. Experiments in Fluids. Volume 25, Issue 3, pp 177–189. <a href="https://doi.org/10.1007/s003480050222">10.1007/s003480050222</a>.<br>
See also: <a href="http://rescience.github.io/ten-years/">Ten Years Reproducibility Challenge<a> - <a href="https://github.com/ReScience/ten-years/issues/1#issuecomment-553313703">paper #24</a>.

Paper: <br>
<a href="https://github.com/quenot/opflow/raw/master/jeif98/jeif98-author.pdf">jeif98-author.pdf</a> : author version of the reference paper. <br>

Data: <br>
<a href="https://github.com/quenot/opflow/raw/master/jeif98/synthetic">synthetic</a> : synthetic images and vector fields for the calibration experiments (23MB). <br>
<a href="https://github.com/quenot/opflow/raw/master/jeif98/real">real</a> : real images for the freezing cavity experiments (1 MB).

Code: <br>
<a href="https://github.com/quenot/opflow/raw/master/jeif98/original">original</a> : source of the original optical flow code. <br>
<a href="https://github.com/quenot/opflow/raw/master/jeif98/threaded">threaded</a> : source of the threaded optical flow code. <br>
There are two versions of the code, "original" which is the version of the code as it was at the time at which the experiments were carried out and "threaded" which differ from the former by the insertion of a single OpenMP parallel directive. Both contain a makefile for compilation on linux systems, either with the GNU or with the intel compiler (if available). It should be easy to modify the top-level makefile to have the code compile on any Unix-like system.

Scripts: <br>
The top-level script for running the "synthetic" experiments is <a href="https://github.com/quenot/opflow/raw/master/jeif98/bat8.sh">bat8.sh</a> </br>
The top-level script for running the "real" experiments is <a href="https://github.com/quenot/opflow/raw/master/jeif98/treal60.sh">treal60.sh</a> </br>
These scripts assume the presence of executables compiled with both the GNU and the Intel compilers. If Intel executables are absent, the corresponding lines in the scripts should be commented out.
