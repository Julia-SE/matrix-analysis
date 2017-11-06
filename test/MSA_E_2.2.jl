include("C:\\Users\\darien.shannon\\Documents\\Code\\Julia\\Julia-SE\\matrix-analysis\\src\\2DTruss.jl")

COORD = [Node(0,0)
         Node(3000/tan(deg2rad(30)) + 3000, 0)
         Node(0, 6000)
         Node(3000/tan(deg2rad(30)) + 3000, 6000)
         Node(3000/tan(deg2rad(30)), 3000)]

MSUP = [Support(1, (1, 1))
        Support(2, (1, 1))
        Support(3, (1, 1))
        Support(4, (1, 1))]

EM = Dict( 1 => Material(200e3))

CP = Dict( 1 => Section(5),
           2 => Section(3))


MPRP = [Element(5, 1, 1, 1)
        Element(5, 3, 1, 1)
        Element(5, 2, 1, 2)
        Element(5, 4, 1, 2)]

memberLoads = []

nodeLoads = [(5, NLoad(1000., 0.))]

forces, reactions, displacements = run_analysis(COORD, MSUP, EM, CP, MPRP, nodeLoads)
println(forces)
println(reactions)
println(displacements)

Base.Test.@test forces ≈ [368.7, 368.7, -255.4, -255.4] atol= .5

Base.Test.@test displacements ≈ [2.55, 0.0] atol = .01
