import Base: endof, setindex!, getindex

immutable PetscVec_s; end
typealias PetscVec_ Ptr{PetscVec_s}

function getindex(x::PetscVec_, ind::Array{Int64})

	vals = zeros( length(ind) );
	ccall( (:VecGetValues, saltlib), Cint, 
		(PetscVec_, Int, Ptr{Int}, Ptr{Cdouble}), 
		x, length(ind), int32(ind-1), vals
	);
	
	if( length(vals) == 1)
		return vals[1];
	else
		return vals;
	end
end

getindex(x::PetscVec_, ind::Integer) = getindex(x, [ind]);
getindex(x::PetscVec_, ind::Range1) = getindex(x, [ind]);

function setindex!(x::PetscVec_, vals::Array{Cdouble}, ind::Array{Int})
    length(vals) != length(ind) && throw(ArgumentError("vals and ind must have same length"))
    
    ccall( (:VecSetValues, saltlib), Cint, 
	  (PetscVec_, Int, Ptr{Int}, Ptr{Cdouble}, Int),
	  x, length(ind), int32(ind-1), vals, 1
	  )  # 1 is INSERT_VALUES
    
    jul
end

setindex!(x::PetscVec_, v::Array{Cdouble}, i::AbstractArray) =
  setindex!(x, v, copy!(Array(Int, length(i)), i))

setindex!(x::PetscVec_, v::AbstractArray, i::AbstractArray) =
  setindex!(x, copy!(Array(Cdouble, length(v)), v), i)

setindex!(x::PetscVec_, v::Real, i::Integer) = 
   setindex!(x, Cdouble[val], Int[i])

function endof(x::PetscVec_)

	N = [0];
	ccall( (:VecGetSize, saltlib), Cint,
		(PetscVec_, Ptr{Int} ), x, N );
	N[1]
end