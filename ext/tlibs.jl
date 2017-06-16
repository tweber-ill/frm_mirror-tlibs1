#
# tlibs Julia module
# @author Tobias Weber <tobias.weber@tum.de>
# @date 23-apr-2017
# @license GPLv2 or GPLv3
#

__precompile__()
module tl

t_real = Float64
enable_debug = 0


#
# initialises tlibs
#
function __init__()
	ccall((:load_tlibs, :tlibs_jl), Void, (Cint,), enable_debug)
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
function fit(fkt, x, y, yerr; fixed = [], values = Dict(), errors = Dict())
	# find number of function arguments
	meth = methods(fkt).ms[1]
	num_args = meth.sig.parameters.length - 1
	num_free_params = num_args - 1


	# get function argument names
	strMeth = repr(meth)
	strArgs = strMeth[searchindex(strMeth, "(")+1 : searchindex(strMeth, ")")-1]
	arrArgs = map(strip, split(strArgs, ","))
	arrArgs = map(String, arrArgs)
	arrParams = arrArgs[2 : length(arrArgs)]	# only free params


	# build values/errors array
	iParam = 0
	arrValues = Array{t_real, 1}(num_free_params)
	arrErrors = Array{t_real, 1}(num_free_params)

	for strArg in arrArgs
		# skip "x" parameter
		if(iParam == 0)
			iParam += 1
			continue
		end

		valArg = get(values, strArg, 0.)
		errArg = get(errors, strArg, valArg*0.1)
		
		arrValues[iParam] = valArg
		arrErrors[iParam] = errArg
		iParam += 1
	end

	# map to a C function pointer with "num_args" arguments
	cfkt = cfunction(fkt, t_real, NTuple{num_args, t_real})

	# call C function pointer
	bOk = ccall((:fit, :tlibs_jl),
		# return type
		Cint,

		# arg types
		(Ptr{Void}, Csize_t,
		Ptr{t_real}, Ptr{t_real}, Ptr{t_real}, Csize_t,
		Array{String, 1}, Array{String, 1},
		Ptr{t_real}, Ptr{t_real}),

		# args
		cfkt, num_free_params, x, y, yerr, length(x), arrArgs, fixed, arrValues, arrErrors)



	# build map of values & errors
	dictRet = Dict()

	for (strParam, valArg, valErr) in zip(arrParams, arrValues, arrErrors)
		dictRet[strParam] = valArg
		dictRet[strParam * "_err"] = valErr
	end

	dictRet["<valid>"] = bOk
	return dictRet
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



fitresult = tl.fit(tl.gauss_model_amp, Es, cts, cts_err, fixed = ["offs"],
	values = Dict("x0" => -1.5, "sigma" => 0.5, "amp" => 100., "offs" => 50.),
	errors = Dict("x0" => 0.25, "sigma" => 0.25, "amp" => 20., "offs" => 10.))
println(fitresult)

# -----------------------------------------------------------------------------
