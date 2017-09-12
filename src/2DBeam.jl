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
    I::Float64
end

struct Element
    NodeB::Int
    NodeE::Int
    Material
    Section
end

struct CLoad
    value::Float64
    distance::Float64
end

struct MLoad
    value::Float64
    distance::Float64
end

struct NLoad{T <: Tuple{Vararg{AbstractFloat}}}
    loads::T
end

NLoad(loads::Vararg{AbstractFloat}) = NLoad(loads)


struct DLoad
    value1::Float64
    value2::Float64
    distance1::Float64
    distance2::Float64
end


struct ALoad
    value::Float64
    distance::Float64
end

struct SurfaceShearLoad
    value::Float64
    distance1::Float64
    distance2::Float64
end



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

function k(Xb, Yb, Xe, Ye, E, I)
    δx = Xe - Xb
    δy = Ye - Yb
    L = sqrt(δx^2 + δy^2)

    (E*I/L^3)*[12. 6L -12 6L
               6L 4L^2 -6L 2L^2
               -12 -6L 12 -6L
               6L 2L^2 -6L 4L^2]
end


function Qf(Xb, Yb, Xe, Ye, cLoad::CLoad)
    δx = Xe - Xb
    δy = Ye - Yb
    L = sqrt(δx^2 + δy^2)
    W = cLoad.value
    l₁ = cLoad.distance
    l₂ = L - cLoad.distance

    FSb = (W*l₂^2/L^3)*(3l₁ + l₂)
    FMb = (W*l₁*l₂^2/L^2)
    FSe = (W*l₁^2/L^3)*(l₁ + 3l₂)
    FMe = (-W*l₁^2*l₂/L^2)

    [FSb; FMb; FSe; FMe]
end




function Qf(Xb, Yb, Xe, Ye, dLoad::DLoad)
    δx = Xe - Xb
    δy = Ye - Yb
    L = sqrt(δx^2 + δy^2)
    w₁ = dLoad.value1
    w₂ = dLoad.value2
    l₁ = dLoad.distance1
    l₂ = L - dLoad.distance2

    isUniform = w₁ == w₂

    FSb, FMb, FSe, FMe, = (0., 0., 0., 0.)

    if isUniform
        w = w₁

        FSb = (w*L/2)*(1 - (l₁/L^4)*(2L^3 - 2l₁^2*L + l₁^3) - (l₂^3/L^4)*(2L - l₂))
        FMb = (w*L^2/12)*(1 - (l₁^2/L^4)*(6L^2 - 8l₁*L + 3l₁^2) - (l₂^3/L^4)*(4L - 3l₂))
        FSe = (w*L/2)*(1 - (l₁^3/L^4)*(2L - l₁) - (l₂/L^4)*(2L^3 - 2l₂^2 * L + l₂^3))
        FMe = (-w*L^2/12)*(1 - (l₁^3/L^4)*(4L - 3l₁) - (l₂^2/L^4)*(6L^2 - 8l₂*L + 3l₂^2))
    else
        FSb = w₁*(L-l₁)^3/(20*L^3)*( (7L + 8l₁) - (l₂*(3L + 2l₁)/(L-l₁))
              *(1 + (l₂/(L-l₁)) + (l₂^2/(L-l₁)^2)) + (2l₂^4/(L-l₁)^3))
              +(w₂*(L-l₁)^3/20L^3)*((3L + 2l₁)*(1 + l₂/(L - l₁)
                + l₂^2/(L-l₁)^2) - (l₂^3/(L-l₁)^2)*(2 + (15L - 8l₂)/(L-l₁)))

        FMb = (w₁*(L-l₁)^3/60L^2)*(3(L+4l₁) - (l₂*(2L + 3l₁)/(L-l₁))*
               (1 + (l₂/(L-l₁)) + (l₂^2/(L-l₁)^2)) + (3l₂^4/(L-l₁)^3))
               + (w₂*(L-l₁)^3/60L^2)*((2L + 3l₁)*(1 + (l₂/(L-l₁))
               + (l₂^2/(L-l₁)^2)) - (3l₂^3/(L-l₁)^2)*(1 + (5L-4l₂)/(L-l₁)))

        FSe = ((w₁ + w₂)/2)*(L - l₁ - l₂) - FSb

        FMe = ((L- l₁ - l₂)/60)*(w₁*(-2L + 2l₁ - l₂) - w₂*(L - l₁ + 2l₂))
               + FSb*L - FMb

    end
    [FSb; FMb; FSe; FMe]
end

function Qf(Xb, Yb, Xe, Ye, loads::Array{CLoad,1})
    sum(map(load -> Qf(Xb, Yb, Xe, Ye, load), loads))
end

function Qf(Xb, Yb, Xe, Ye, loads::Array{DLoad,1})
    sum(map(load -> Qf(Xb, Yb, Xe, Ye, load), loads))
end

function Qf(Xb, Yb, Xe, Ye, loads::Array{Any,1})
    if size(loads, 1) == 0
        return [0.0; 0.0; 0.0; 0.0]
    else
        sum(map(load -> Qf(Xb, Yb, Xe, Ye, load), loads))
    end
end

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
            if MSUP[supportIndex].Releases[1] == 1             #     if node is restrained in the X direction...
                 NSC[2*nodeIndex - 1] = restraintCoord           #       NSC X value for this node is the next restraint coordinate
                restraintCoord += 1
            else                                                #     otherwise, this node is free in the X direction
                 NSC[2*nodeIndex - 1] = freeCoord         #       NSC X value for this node is the next free coordinate
                freeCoord += 1
            end
            if MSUP[supportIndex].Releases[2] == 1             #     if node is restrained in the Y direction...
                NSC[2*nodeIndex] = restraintCoord         #       NSC Y value for this node is the next restraint coordinate
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
        I = CP[memberInfo.Section].I
        GK = k(Xb, Yb, Xe, Ye, E, I)

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




function p_pf(COORD, MSUP, NSC, memberLoads, nodeLoads, NCJT)
    NJ = size(COORD, 1)
    NR = sum(map(s -> sum(s.Releases), MSUP))
    numDOF = NJ * NCJT - NR
    NML = length(memberLoads)
    NJL = length(nodeLoads)

    P = zeros(numDOF)::Array{Float64,1}
    Pf = zeros(numDOF)::Array{Float64,1}

    dofIndecies = find(x -> x <= numDOF, NSC)

    for i in 1:size(dofIndecies,1)

        dofIndex = dofIndecies[i]
        remainder = dofIndex%NCJT
        dofOfNode = (remainder == 0 ? NCJT : remainder)
        node = Int(1 + (dofIndex - dofOfNode)/NCJT)


        _nodeLoads =  map(load -> load[2], filter(load -> load[1] == node, nodeLoads))

         for load in _nodeLoads
            P[i] = P[i] + load.loads[dofOfNode]
         end


        for j in 1:size(MPRP,1)
            member = MPRP[j]
            nodeBegin, nodeEnd = member.NodeB, member.NodeE
            Xb = COORD[nodeBegin].X
            Yb = COORD[nodeBegin].Y
            Xe = COORD[nodeEnd].X
            Ye = COORD[nodeEnd].Y

            _memberLoads = map(load -> load[2], filter(load -> load[1] == j, memberLoads))
            _Qf = Qf(Xb, Yb, Xe, Ye, _memberLoads)

            if nodeBegin == node
                Pf[i] = Pf[i] + _Qf[dofOfNode]
            end
            if nodeEnd == node
                Pf[i] = Pf[i] + _Qf[2 + dofOfNode]
            end
        end
    end
    P, Pf
end



function get_d(S, P, Pf)
    inv(S)*(P - Pf)
end



function forces_reactions(COORD, NSC, MSUP, EM, CP, MPRP, d, NCJT)
    NJ = size(COORD, 1)
    NR = sum(map(s -> sum(s.Releases), MSUP))
    numDOF = NJ * NCJT - NR

    R = zeros(NR)
    forces = Array{Array{Float64,1}}(size(MPRP, 1))

    for member in 1:size(MPRP,1)
        memberInfo = MPRP[member]
        nodeBegin, nodeEnd = memberInfo.NodeB, memberInfo.NodeE
        Xb = COORD[nodeBegin].X
        Yb = COORD[nodeBegin].Y
        Xe = COORD[nodeEnd].X
        Ye = COORD[nodeEnd].Y
        E = EM[memberInfo.Material].E
        I = CP[memberInfo.Section].I

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

        BK = k(Xb, Yb, Xe, Ye, I, E)

        _memberLoads = map(load -> load[2], filter(load -> load[1] == member, memberLoads))
        _Qf = Qf(Xb, Yb, Xe, Ye, _memberLoads)

        Q = BK * U + _Qf
        forces[member] = Q

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
    forces, R
end


function run_analysis(COORD, MSUP, EM, CP, MPRP, memberLoads, nodeLoads)
    NCJT = 2
    NSC = nsc(COORD, MSUP, NCJT)
    GS = gs(COORD, MSUP, MPRP, EM, CP, NSC, NCJT)
    P, Pf = p_pf(COORD, MSUP, NSC, memberLoads, nodeLoads, NCJT)
    d = get_d(GS, P, Pf)
    forces, reactions = forces_reactions(COORD, NSC, MSUP, EM, CP, MPRP, d, NCJT)
end
