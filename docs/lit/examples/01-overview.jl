#=
# [MRIPulses overview](@id 01-overview)

This page illustrates the Julia package
[`MRIPulses`](https://github.com/MagneticResonanceImaging/MRIPulses.jl).
=#

#srcURL


# ## Setup

# Packages needed here.

using MRIPulses
using MIRTjim: jim, prompt
using InteractiveUtils: versioninfo


# The following line is helpful when running this file as a script;
# this way it will prompt user to hit a key after each figure is displayed.

isinteractive() ? jim(:prompt, true) : prompt(:draw);


# ## Overview

#=
See the other demo(s).
=#

include("../../../inc/reproduce.jl")
