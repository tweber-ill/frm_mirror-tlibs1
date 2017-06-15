#
# tlibs Julia module
# @author Tobias Weber <tobias.weber@tum.de>
# @date 23-apr-2017
# @license GPLv2 or GPLv3
#

__precompile__()
module tl

t_real = Float64


#
# initialises tlibs
#
function __init__()
	ccall((:load_tlibs, :tlibs_jl), Void, ())
end


#
# loads instrument data files
#
function loadinstr(strFile::String) :: Array{Any, 1}
	scandata = ccall((:load_instr, :tlibs_jl), Array{Any, 1}, (Cstring,), strFile)

	thedict = Dict{String,String}()
	for (key, val) in zip(scandata[3], scandata[4])
		thedict[key] = val
	end

	# header names, data matrix, property dict
	return [ scandata[1], scandata[2], thedict ]
end



#
# fitting
#
function fit(fkt, x, y, yerr)
	# find number of function arguments
	meth = methods(fkt).ms[1]
	num_args = meth.sig.parameters.length - 1
	num_free_params = num_args - 1

	# map to a C function pointer
	if(num_args == 2)
		cfkt = cfunction(fkt, t_real, (t_real, t_real))
	elseif (num_args == 3)
		cfkt = cfunction(fkt, t_real, (t_real, t_real, t_real))
	elseif (num_args == 4)
		cfkt = cfunction(fkt, t_real, (t_real, t_real, t_real, t_real))
	elseif (num_args == 5)
		cfkt = cfunction(fkt, t_real, (t_real, t_real, t_real, t_real, t_real))
	elseif (num_args == 6)
		cfkt = cfunction(fkt, t_real, (t_real, t_real, t_real, t_real, t_real, t_real))
	elseif (num_args == 7)
		cfkt = cfunction(fkt, t_real, (t_real, t_real, t_real, t_real, t_real, t_real, t_real))
	elseif (num_args == 8)
		cfkt = cfunction(fkt, t_real, (t_real, t_real, t_real, t_real, t_real, t_real, t_real, t_real))
	elseif (num_args == 9)
		cfkt = cfunction(fkt, t_real, (t_real, t_real, t_real, t_real, t_real, t_real, t_real, t_real, t_real))
	else
		println("Invalid or unsupported number of arguments for fit function.")
		return false
	end

	# call C function pointer
	bOk = ccall((:fit, :tlibs_jl),
		Cint,	# return type
		(Ptr{Void}, Csize_t,
		Ptr{t_real}, Ptr{t_real}, Ptr{t_real}, Csize_t),	# arg types
		cfkt, num_free_params, x, y, yerr, length(x))

	return bOk != 0
end


#
# Gauss model
#
function gauss_model_amp(x, x0, sigma, amp, offs)
	return amp * exp(-0.5 * ((x-x0)/sigma)*((x-x0)/sigma)) + offs
end

#
# Lorentz model
#
function lorentz_model_amp(x, x0, hwhm, amp, offs)
	return amp*hwhm*hwhm / ((x-x0)*(x-x0) + hwhm*hwhm) + offs
end



end 	# tl




# -----------------------------------------------------------------------------
# test


a = tl.loadinstr("/home/tweber/Measurements/mira-mgv2o4-17/data/11797_00025700.dat")
cols = a[1]
data = a[2]
#println(cols)
#println(typeof(data))

idxE = findfirst(cols, "E")
idxCtr = findfirst(cols, "ctr1")

Es = data[:,idxE]
cts = data[:,idxCtr]
cts_err = sqrt(cts)

println(lpad(cols[idxE], 12), " ", lpad(cols[idxCtr], 12))
for (E, ct) in zip(Es, cts)
	println(lpad(E, 12), " ", lpad(ct, 12))
end



tl.fit(tl.gauss_model_amp, Es, cts, cts_err)


# -----------------------------------------------------------------------------
