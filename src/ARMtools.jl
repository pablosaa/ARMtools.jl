"""
Julia Package including functions to read and process data from the DOE ARM program.

The package is mainly focused to work with instruments from the NSA facility in Utqiagvik, Alaska.
Nonetheless since most of the ARM instrumentation is standard, the package can also be used for
other facilities.

(c) 2020, Pablo Saavedra Garfias
University of Leipzig
Faculty of Physics and Geosciences

See LICENSE
"""
module ARMtools

using NCDatasets
# * Read INTERPOLATE radiosonde tools:
"""
Function getSondeData(file_name::String)

This will read the ARM datafile 'file_name.nc' and return a Dictionary
with the defaul data fields.
The default data fields are:
* :time
* :heitht
* :T
* :Pa
* :RH
* :qv
* :U
* :V
* :WD
* :WS
* :θ

"""
function getSondeData(sonde_file::String; addvars=(), onlyvars=() )

    ncvars = Dict(:time=>"time",
                  :height=>"height",
                  :T=>"temp",
                  :Pa=>"bar_pres",
                  :RH=>"rh",
                  :qv=>"sh",
                  :U=>"u_wind",
                  :V=>"v_wind",
                  :WSPD=>"wspd",
                  :WDIR=>"wdir",
                  :θ=>"potential_temp")
    # for INTERPOLATE RS data:
    if ~isempty(addvars)
        # finding whether the variables are already as default:
        idx = map(x->findfirst(==(x), ncvars), addvars) .|> isnothing |> findall
        # converting to Dict the added variables:
        tmp = map(x->Symbol(uppercase(x))=>x, addvars[idx]) |> Dict
        # merging with default variables:
        ncvars = merge(ncvars, tmp)
    end
    if ~isempty(onlyvars)
        # converting only variables to Dict:
        tmp = map(x->Symbol(uppercase(x))=>x, onlyvars[idx]) |> Dict
        ncvars = merge(Dict(:time=>"time", :height=>"height"), ncvars[3:end])
    end

    output = Dict()
    isfile(sonde_file) ? "reading $sonde_file" : "file not found!"
    ncin = NCDataset(sonde_file)
    for var ∈ ncvars
        str_var = var[2]
        println(str_var)
        tmp_var = ncin[str_var][:,:]
        if haskey(ncin[str_var].attrib, "missing_value")
            miss_val = ncin[str_var].attrib["missing_value"]
            tmp_var[tmp_var .≈ miss_val] .= NaN
        end
        # filling the output variable:
        key_var = var[1]
        output[key_var] = tmp_var
    end
    close(ncin)
    return output
    
    ##T_C  = copy(ncin["temp"][:,:])
    ##rh = copy(ncin["rh"][:,:])
    ##time = ncin["time"][:]
    ##height= copy(ncin["height"][:]) # km
    ##wind = copy(ncin["wspd"][:,:])  # m/s
    ##miss_val = ncin["wspd"].attrib["missing_value"]
    ##wind[wind .≈ miss_val] .= NaN
    ##wdir = copy(ncin["wdir"][:,:])  # deg
    ##miss_val = ncin["wdir"].attrib["missing_value"]
    ##wdir[wdir .≈ miss_val] .= NaN
    ##
    ##uwind = ncin["u_wind"][:,:]; # m/s
    ##miss_val = ncin["u_wind"].attrib["missing_value"]
    ##uwind[uwind .≈ miss_val] .= NaN
    ##
    ##vwind = ncin["v_wind"][:,:]; # m/s
    ##miss_val = ncin["v_wind"].attrib["missing_value"]
    ##vwind[vwind .≈ miss_val] .= NaN
    ##
    #### Potential temperature
    ##θ = ncin["potential_temp"][:,:] # K
    ##miss_val = ncin["potential_temp"].attrib["missing_value"]
    ##θ[θ .≈ miss_val] .= NaN
    ##
    #### New variables for IVT
    ##qv = ncin["sh"][:,:];  # gr/gr
    ##miss_val = ncin["sh"].attrib["missing_value"]
    ##qv[qv .≈ miss_val] .= NaN
    ##
    ##pa = ncin["bar_pres"][:,:];  # kPa
    ##miss_val = ncin["bar_pres"].attrib["missing_value"]
    ##pa[pa .≈ miss_val] .= NaN
    ##
    ##close(ncin)
    ##return Dict(:time=>time, :height=>height[:,1], :T=>T_C) #, U=uwind, V=vwind, RH=rh, WS=wind, WD=wdir, QV=qv, Pa=pa, θ = θ)

end


end # module
# Main file containing the package module
# See LICENSE
