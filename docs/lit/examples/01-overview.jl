#=
# [MRIPulses overview](@id 01-overview)

This page illustrates the Julia package
[`Template`](https://github.com/MagneticResonanceImaging/MRIPulses.jl).

This page was generated from a single Julia file:
[01-overview.jl](@__REPO_ROOT_URL__/01-overview.jl).
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
=#

include("../../../inc/reproduce.jl")
