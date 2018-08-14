# Fast track simulation extension for the Acts project

This project provides the code to run a fast track simulation on top of 
`Acts::Core`.

Being a fast track simulation, it uses the reconstruction geometry description,
i.e. the `Acts::TrackingGeometry` as a simulation geometry with simplified and
particle parameterised material effects.

Dependencies for the Core components are:
    * `Acts::Core` and consequently `Eigen` and `Boost`

Optional dependency exists for 
    * `Geant4` for optional functionality taken from the full simulation toolkit
    

