immutable PetscVec_s; end
typealias PetscVec_ Ptr{PetscVec_s}

function getindex(x::PetscVec_, ind)

	if typeof(ind) <: Integer || typeof(ind) <: Range1
		ind = [ind]
	end

	vals = zeros( length(ind) );
	ccall( (:VecGetValues, saltlib), Cint, 
		(PetscVec_, Cint, Ptr{Cint}, Ptr{Cdouble}), 
		x, length(ind), int32(ind-1), vals
	);
	
	if( length(vals) == 1)
		return vals[1];
	else
		return vals;
	end
end

function setindex!(x::PetscVec_, vals, ind)

	length(vals) != length(ind) && throw(ArgumentError("vals and ind must have same length"))

	if typeof(ind) <: Integer || typeof(ind) <: Range1{Int64}
		ind = [ind]
	end

	if typeof(vals) <: Cdouble
		vals = [vals]
	end

	ccall( (:VecSetValues, saltlib), Cint, 
		(PetscVec_, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint),
		x, length(ind), int32(ind-1), vals, 1
	);  # 1 is INSERT_VALUES
end

import Base.endof  # for some reason you don't have to do this for getindex and setindex

function endof(x::PetscVec_)

	N = [0];
	ccall( (:VecGetSize, saltlib), Cint,
		(PetscVec_, Ptr{Cint} ), x, N );
	N[1]
end