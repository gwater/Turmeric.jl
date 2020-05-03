# This file is MIT licensed

module IntervalRootFinding2

using ForwardDiff
using StaticArrays

using LinearAlgebra: I, Diagonal

import Base: ⊆, show, big, \

include("helpers.jl")
include("contractors.jl")

const default_tolerance = 1e-7
const default_contractor = Krawczyk()

include("parallel.jl")
include("roots.jl")

end
