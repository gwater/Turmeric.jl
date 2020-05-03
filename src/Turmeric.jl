module Turmeric

include("helpers.jl")
include("contractors.jl")

const default_tolerance = 1e-7
const default_contractor = Krawczyk()

include("threadsearch.jl")
include("roots.jl")
include("contract.jl")

end # module
