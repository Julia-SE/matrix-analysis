include("C:\\Users\\darien.shannon\\Documents\\Code\\Julia\\JuliaSE\\src\\2DBeam.jl")

COORD = [Node(0, 0)
         Node(6, 0)
         Node(10, 0)
         Node(20, 0)]

MSUP = [Support(1, 1, 1)
        Support(2, 0, 0)
        Support(3, 1, 0)
        Support(4, 1, 0)]

EM = Dict( 1 => Material(28e6))

CP = Dict( 1 => Section(5.8e9/(1000^4)),
           2 => Section(1.5*5.8e9/(1000^4)))


MPRP = [Element(1, 2, 1, 2)
        Element(2, 3, 1, 1)
        Element(3, 4, 1, 1)]

memberLoads = [(1, DLoad(30, 0, 0, 6))
               (3, CLoad(150, 5))]

nodeLoads = [(2, NLoad((-200, 0)))
             (3, NLoad((0, -90)))]



forces, reactions = run_analysis(COORD, MSUP, EM, CP, MPRP, memberLoads, nodeLoads)
println(forces)
println(reactions)
