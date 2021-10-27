using DelimitedFiles

"""
ExtendedUltrasphericalArc{T}()

This extends the UltrasphericalArc(0) domain as viewed in Cartesian coordinates. The original
domain is [-1,1]. My mapping the function to the sum space 
    
    ExtendedChebyshevT + Weighted(ChebyshevU())

and mapping back for values x ̸∈ [-1,1], we can define natural extensions of the arc polynomials. 


"""

struct ExtendedUltrasphericalArc{T, RR} <: Basis{T} 
    R::RR
end
ExtendedUltrasphericalArc{T}(R::RR) where {T,RR} = ExtendedUltrasphericalArc{T,RR}(R)

function ExtendedUltrasphericalArc{T}() where T
    R = Matrix{T}(undef,31,31)
    ExtendedUltrasphericalArc{T}(R)
end

ExtendedUltrasphericalArc() = ExtendedUltrasphericalArc{Float64}()

axes(P::ExtendedUltrasphericalArc) = (Inclusion(ℝ), _BlockedUnitRange(1:2:∞))

==(a::ExtendedUltrasphericalArc, b::ExtendedUltrasphericalArc) = true

function getindex(P::ExtendedUltrasphericalArc{T}, x::Real, j::Int)::T where T
    U = UltrasphericalArc()
    x in ChebyshevInterval() && return U[CircleCoordinate(acos(x)), j]
    
    R = readdlm("R.txt")

    c_p = Array{Float64}(undef,size(R,1))
    c_p .= 0
    c_p[j] = 1
    
    # Coefficient conversion to sum space
    c = R \ c_p
    eT = ExtendedChebyshevT()
    return sum(c[1:2:end].*eT[x,1:(size(R,1)+1)÷2])

    # if P.f == 0
    #     xy = axes(U,1)
    #     x = first.(xy)
    #     P.R[:,:].=0
    #     P.R[:,1] = (U \ broadcast(x -> ChebyshevT()[x,1], x))[1:size(P.R,1)]
    #     @showprogress for j = 1:size(P.R,2)÷2
    #         P.R[:,2j+1] = (U \ broadcast(x -> ExtendedChebyshevT()[x,j+1], x))[1:size(P.R,1)]
    #         P.R[:,2j] = (U \ broadcast(x -> sqrt(1-x^2)*ChebyshevU()[x,j], x))[1:size(P.R,1)]
    #     end
    #     P.f = 1
    # end
end
