function plotTEslice(y)

	if !isdefined(:plot)
		using PyPlot;
	end

	plot(y);

end