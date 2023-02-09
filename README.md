# MuscRat.jl

A Muon Simulation for Cosmic Ray Analysis Tasks

<!-- Couldn't figure out how to get markdown to resize the drawing w/o going to full HTML -->
<!-- ![A cute muskrat](docs/640px-Bisamratte-drawing.jpg) -->
<p align="left">
    <img src="docs/640px-Bisamratte-drawing.jpg", width="160" title="A cute muskrat">
</p>

The MuscRat project is designed to analyze the flux of cosmic ray muons in detectors.
This code runs in Julia 1.x (tested for 1.7 and up).

Analytic geometry allows MuscRat to compute the flux through muon detectors of simple, regular shape. Shapes implmented are:
* a horizontal cylinder (`Hcylinder`)
* a rectangular prism (`Box`)

The muon flux is generated based on one of the two analytic formulas found in Su et al (2021): those of Reyna 2006 or Chatzidakis et al. 2015. At the moment, these are only sea-level values.

## References

* Su, N., et al., (2021). "A Comparison of Muon Flux Models at Sea Level for Muon Imaging and Low Background Experiments." _Frontiers in Energy Research_, **9**. [doi:10.3389/fenrg.2021.750159](https://doi.org/10.3389/fenrg.2021.750159)
* Reyna, D. (2006). "A Simple Parameterization of the Cosmic-Ray Muon Momentum Spectra at the Surface as a Function of Zenith Angle" [arXiv:hep-ph/0604145](https://arxiv.org/abs/hep-ph/0604145)
* Chatzidakis, S., Chrysikopoulou, S., Tsoukalas, L. H. (2015). "Developing a cosmic ray muon sampling capability for muon tomography and monitoring applications." _Nuclear Instruments and Methods in Physics Research, Section A_, **804**. [doi:10.1016/j.nima.2015.09.033](https://doi.org/10.1016/j.nima.2015.09.033)