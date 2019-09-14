# porousTransportUpscalingFoam
OpenFOAM solver for the formal upscaling of transport with surface reactions in porous media.

[![DOI](https://zenodo.org/badge/208470601.svg)](https://zenodo.org/badge/latestdoi/208470601)

Requirements
------------

_OpenFOAM-6_ and __Matlab_.

Usage
-----

Simply run the _Allwmake_ script to compile the solver and the library with additional boundary conditions.

Content
-------

* _transportUpscalingFoam_: solver to calculate effective transport coefficients from a reactive transport equation. This effectively solves the cell problems (see references).

* _libUpscalingBCs_: library containing the required new boundary conditions for _transportUpscalingFoam_.

* _tutorials_: contains the _packingCell_ tutorial (main example application) and results from the work cited in the references.

* _etc_ : data and _Matlab_ scripts (_Chebfun_) for verification.

Authors
-------

[Federico Municchi](https://github.com/fmuni) 

[Matteo Icardi](https://github.com/matteoicardi)

References
----------

[Municchi, Federico, and Matteo Icardi. "Macroscopic models for heterogeneous reactions in porous media." arXiv preprint arXiv:1909.02818 (2019)](https://arxiv.org/abs/1909.02818)

Aknowledgements
---------------

This work has been funded by the European Union's Horizon 2020
research and innovation programme, grant agreement number 764531, "SECURe -- Subsurface Evaluation of Carbon capture and storage and Unconventional risks".
