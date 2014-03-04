immutable Mode_s; end
typealias Mode_ Ptr{Mode_s}
DestroyMode(m::Mode_) = ccall((:DestroyMode, saltlib), Void, (Mode_,), m )

type Mode
	m::Mode_
	psi::PetscVec_
	function Mode(m::Mode_)
		md = new(
			m,
			ccall( (:GetVpsi, saltlib), PetscVec_, (Mode_,), m)
		);
		finalizer(md, DestroyMode)
		return md
	end

end