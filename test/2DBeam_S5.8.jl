include("C:\\Users\\darien.shannon\\Documents\\Code\\Julia\\JuliaSE\\src\\2DBeam.jl")

COORD = [Node(0,0)
         Node(120, 0)
         Node(360, 0)
         Node(480, 0)]

MSUP = [Support(1, 1, 0)
        Support(2, 1, 0)
        Support(3, 1, 0)
        Support(4, 1, 1)]

EM = Dict( 1 => Material(29e3))

CP = Dict( 1 => Section(350),
           2 => Section(500))


MPRP = [Element(1, 2, 1, 1)
        Element(2, 3, 1, 1)
        Element(3, 4, 1, 2)]

memberLoads = [(2, DLoad(0.1667, 0.1667, 0, 120))
               (2, CLoad(25, 180))
               (3, DLoad(0.25, 0, 0, 120))]

nodeLoads = [(1, NLoad((0, -480)))]



forces, reactions = run_analysis(COORD, MSUP, EM, CP, MPRP, memberLoads, nodeLoads)
println(forces)
println(reactions)
