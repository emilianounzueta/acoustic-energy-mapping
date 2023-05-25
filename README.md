# Acoustic Energy Mapping

This repository cointains the code for processing acoustic signals captured by microphones for generating an acoustic energy map in a 3D-plot and as a heatmap with simulated and real signals.

## Table of Contents
- [Description](#description)
- [Papers](#papers)
- [License](#license)

## Description
This repository cointains the code that can be used for the generation of an acoustic energy map by processing signals captured by an acoustic sensor network. The concepts are explained in the article [On the challenges of Acoustic Energy Mapping Using a WASN: Synchronization and Audio Capture](https://www.mdpi.com/1424-8220/23/10/4645). For the processing of the signals it is required to know: the position of the microphones, the acoustic signal emmited by the source and the position of the "proof points" where the energy will be calculated. 

![Image 1](/heat1.png)
![Image 2](/map.png)

The repository is divided in three sections: 
1. Characterization. Where the signals are simulated and a number of variables in the setup of the WASN, the "proof points" and the input signal of the source can be modified. 
2. Recording processing. For the processing of real acoustic signals captured by a 4-node WASN, this includes the post-processing of the signals and the implementation of the acoustic energy mapping model.
3. Recordings. A set of recordings for testing.

## Papers

This code is result of the article published in "Sensors"
journal.

> García-Unzueta, E.E.; Mendez-Monroy, P.E.; Rascon, C. 
> On the Challenges of Acoustic Energy Mapping Using a 
> WASN: Synchronization and Audio Capture. Sensors 2023, 23, 4645
> https://doi.org/10.3390/s23104645

## License

See license.txt in the root directory of the repository.
