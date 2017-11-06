struct Node
    X::Float64
    Y::Float64
end

struct Support{T <: Tuple{Vararg{Int}}}
    Node::Int
    Releases::T
end

struct Material
    E::Float64
end

struct Section
    Area::Float64
end

struct Element
    NodeB::Int
    NodeE::Int
    Material
    Section
end


struct NLoad{T <: Tuple{Vararg{AbstractFloat}}}
    loads::T
end

NLoad(loads::Vararg{AbstractFloat}) = NLoad(loads)

function gK(Xb, Yb, Xe, Ye, E, A)
  δx = Xe - Xb
  δy = Ye - Yb
  L = sqrt(δx^2 + δy^2)
  _cos = δx/L
  _sin = δy/L
  _cos² = _cos * _cos
  _sin² = _sin * _sin
  _cos_sin = _cos * _sin
  (A*E/L)*[_cos² _cos_sin -_cos² -_cos_sin
   _cos_sin _sin² -_cos_sin -_sin²
   -_cos² -_cos_sin _cos² _cos_sin
   -_cos_sin -_sin² _cos_sin _sin²]
end

function T(Xb, Yb, Xe, Ye)
    δx = Xe - Xb
    δy = Ye - Yb
    L = sqrt(δx^2 + δy^2)
    _cos = δx/L
    _sin = δy/L

    [_cos _sin 0 0
     -_sin _cos 0 0
     0 0 _cos _sin
     0 0 -_sin _cos]
end

function k(Xb, Yb, Xe, Ye, E, A)
    δx = Xe - Xb
    δy = Ye - Yb
    L = sqrt(δx^2 + δy^2)

    (E*A/L)*[1. 0 -1 0
            0 0 0 0
            -1 0 1 0
            0 0 0 0]
end


COORD = [Node(0, 0)
         Node(288, 0)
         Node(576, 0)
         Node(864, 0)
         Node(288, 216)
         Node(576, 216)]

MSUP = [Support(1, (1, 1))
        Support(3, (0, 1))
        Support(4, (0, 1))]

EM = Dict( 1 => Material(29000),
           2 => Material(10000))

CP = Dict( 1 => Section(8),
           2 => Section(12),
           3 => Section(16))


MPRP = [Element(1, 2, 1, 1)
        Element(2, 3, 1, 1)
        Element(3, 4, 2, 3)
        Element(5, 6, 1, 1)
        Element(2, 5, 1, 1)
        Element(3, 6, 1, 1)
        Element(1, 5, 1, 2)
        Element(2, 6, 1, 2)
        Element(3, 5, 1, 2)
        Element(4, 6, 2, 3)]

loads = [(2, NLoad(0.0, -75.))
         (5, NLoad(25., 0.))
         (6, NLoad(0., -60.))]

#=
COORD = [Node(12*12, 16*12)
         Node(0, 0)
         Node(12*12, 0)
         Node(24*12, 0)]

MSUP = [Support(2, 1, 1)
        Support(3, 1, 1)
        Support(4, 1, 1)]

EM = Dict( 1 => Material(29000))

CP = Dict( 1 => Section(8),
           2 => Section(6))

MPRP = [Element(2, 1, 1, 1)
        Element(3, 1, 1, 2)
        Element(4, 1, 1, 1)]

loads = Dict(1 => (150, -300))
=#
#=
COORD = [Node(0,0)
         Node(10,0)
         Node(0,8)
         Node(6,8)]

MSUP = [Support(1, 1, 1)
        Support(2, 1, 1)
        Support(3, 1, 0)]

EM = Dict( 1 => Material(70e6))

CP = Dict( 1 => Section(.004))

MPRP = [Element(1, 3, 1, 1)
        Element(3, 4, 1, 1)
        Element(1, 4, 1, 1)
        Element(2, 3, 1, 1)
        Element(2, 4, 1, 1)]

loads = Dict(3 => (0.0, -400.0),
             4 => (800.0, -400.0))
=#


function nsc(COORD, MSUP, NCJT)
    #Determine number of joints (NJ)
    NJ = size(COORD, 1)

    #Determine Number of Restraints (NR)
    NR = sum(map(s -> sum(s.Releases), MSUP))

    numDOF = NJ * NCJT - NR

    NSC = Array{Int}(2 * NJ)                                  # Create the NSC with uninitialized values

    freeCoord = 1
    restraintCoord = numDOF + 1
    for nodeIndex in 1:size(COORD, 1)                           # Loop through all of the COORD
        supportIndexes = [s.Node for s in MSUP]
        if in(nodeIndex, supportIndexes)                        #   if this node is a support then...
            supportIndex = findfirst(supportIndexes, nodeIndex) #     get the support Index for this node
            if MSUP[supportIndex].Releases[1] == 1                 #     if node is restrained in the X direction...
                 NSC[2*nodeIndex - 1] = restraintCoord          #       NSC X value for this node is the next restraint coordinate
                restraintCoord += 1
            else                                                #     otherwise, this node is free in the X direction
                 NSC[2*nodeIndex - 1] = freeCoord               #       NSC X value for this node is the next free coordinate
                freeCoord += 1
            end
            if MSUP[supportIndex].Releases[2] == 1                 #     if node is restrained in the Y direction...
                NSC[2*nodeIndex] = restraintCoord               #       NSC Y value for this node is the next restraint coordinate
                restraintCoord += 1
            else                                                #     otherwise, this node is free in the Y direction
                 NSC[2*nodeIndex] = freeCoord             #       NSC Y value for this node is the next free coordinate
                freeCoord += 1
            end
        else                                                    #   otherwise, this node is free
             NSC[2*nodeIndex - 1] = freeCoord             #     NSC X value for this node is the next free coordinate
             NSC[2*nodeIndex] = freeCoord + 1             #     NSC Y value for this node is the next free coordinate
            freeCoord += 2
        end
    end
    NSC
end

#NCJT = 2
#NSC = nsc(COORD, MSUP, NCJT)

function gs(COORD, MSUP, MPRP, EM, CP, NSC, NCJT)
    #Determine number of joints (NJ)
    NJ = size(COORD, 1)

    #Determine Number of Restraints (NR)
    NR = sum(map(s -> sum(s.Releases), MSUP))

    numDOF = NJ * NCJT - NR

    S = zeros(numDOF, numDOF)

    for memberIndex in 1:size(MPRP, 1)
        memberInfo = MPRP[memberIndex]
        nodeBegin, nodeEnd = memberInfo.NodeB, memberInfo.NodeE
        Xb = COORD[nodeBegin].X
        Yb = COORD[nodeBegin].Y
        Xe = COORD[nodeEnd].X
        Ye = COORD[nodeEnd].Y
        E = EM[memberInfo.Material].E
        A = CP[memberInfo.Section].Area
        GK = gK(Xb, Yb, Xe, Ye, A, E)

        NSC_Indecies =
        vcat([NSC[(nodeBegin - 1)*NCJT + x] for x in 1:NCJT],
             [NSC[(nodeEnd - 1)*NCJT + x] for x in 1:NCJT])

        iCount = 0
        for i in NSC_Indecies
            if i == NSC_Indecies[1]
                iCount = 0
            end
            iCount += 1
            if i <= numDOF
                jCount = 0
                for j in NSC_Indecies
                    jCount +=1
                    if j <= numDOF
                         S[j,i] = S[j,i] + GK[jCount, iCount]
                    end
                end
            end
        end
    end
    S
end

function p(COORD, MSUP, NSC, loads, NCJT)
    NJ = size(COORD, 1)
    NR = sum(map(s -> sum(s.Releases), MSUP))
    numDOF = NJ * NCJT - NR
    NJL = length(loads)

    P = zeros(numDOF)::Array{Float64,1}

    dofIndecies = find(x -> x <= numDOF, NSC)

    for i in 1:size(dofIndecies,1)

        dofIndex = dofIndecies[i]
        remainder = dofIndex%NCJT
        dofOfNode = (remainder == 0 ? NCJT : remainder)
        node = Int(1 + (dofIndex - dofOfNode)/NCJT)


        _nodeLoads =  map(load -> load[2], filter(load -> load[1] == node, loads))

         for load in _nodeLoads
            P[i] = P[i] + load.loads[dofOfNode]
         end
    end
    P
end

function get_d(S, P)
    inv(S)*P
end


function forces_reactions_displacements(COORD, NSC, MSUP, EM, CP, MPRP, d, NCJT)
    NJ = size(COORD, 1)
    NR = sum(map(s -> sum(s.Releases), MSUP))
    numDOF = NJ * NCJT - NR

    R = zeros(NR)
    forces = Array{Float64}(size(MPRP,1))

    for member in 1:size(MPRP,1)
        memberInfo = MPRP[member]
        nodeBegin, nodeEnd = memberInfo.NodeB, memberInfo.NodeE
        Xb = COORD[nodeBegin].X
        Yb = COORD[nodeBegin].Y
        Xe = COORD[nodeEnd].X
        Ye = COORD[nodeEnd].Y
        E = EM[memberInfo.Material].E
        A = CP[memberInfo.Section].Area

        v = zeros(NCJT*2)

        NSC_Indecies =
        vcat([NSC[(nodeBegin - 1)*NCJT + x] for x in 1:NCJT],
             [NSC[(nodeEnd - 1)*NCJT + x] for x in 1:NCJT])

        for i in 1:(NCJT*2)
            NSC_Index = NSC_Indecies[i]
            if NSC_Index <= numDOF
                v[i] = d[NSC_Index]
            end
        end

        _T = T(Xb, Yb, Xe, Ye)

        U = _T * v

        BK = k(Xb, Yb, Xe, Ye, E, A)

        Q = BK * U
        forces[member] = Q[3]

        F = _T'Q

        for i in 1:(NCJT*2)
            nscIndex =
                if i<=NCJT
                    (nodeBegin -1)*NCJT + i
                else
                    (nodeEnd - 1)*NCJT + (i - NCJT)
                end
            N = NSC[nscIndex]
            if N > numDOF
                R[N - numDOF] = R[N - numDOF] + F[i]
            end
        end
    end
    forces, R, d
end

function run_analysis(COORD, MSUP, EM, CP, MPRP, loads)
    NCJT = 2
    NSC = nsc(COORD, MSUP, NCJT)
    GS = gs(COORD, MSUP, MPRP, EM, CP, NSC, NCJT)
    P = p(COORD, MSUP, NSC, loads, NCJT)
    d = get_d(GS, P)
    forces, reactions, displacment = forces_reactions_displacements(COORD, NSC, MSUP, EM, CP, MPRP, d, NCJT)
end
