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

The muon flux is generated based on the analytic formulas of ...
