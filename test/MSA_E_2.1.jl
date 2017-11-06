include("C:\\Users\\darien.shannon\\Documents\\Code\\Julia\\Julia-SE\\matrix-analysis\\src\\2DTruss.jl")

COORD = [Node(0,0)
         Node(9000,0)
         Node(6000,4000)]

MSUP = [Support(1, (1, 1))
        Support(2, (1, 1))]

EM = Dict( 1 => Material(200e3))

CP = Dict( 1 => Section(6),
           2 => Section(8))


MPRP = [Element(1, 3, 1, 1)
        Element(2, 3, 1, 2)]

memberLoads = []

nodeLoads = [(3, NLoad(500., 0.))]



forces, reactions, displacements = run_analysis(COORD, MSUP, EM, CP, MPRP, nodeLoads)
println(forces)
println(reactions)
println(displacements)

Base.Test.@test displacements ≈ [2.41, 0.72] atol=.01
