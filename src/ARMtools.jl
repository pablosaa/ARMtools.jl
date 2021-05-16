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

using NCDatasets, Statistics

## * Auxiliary functions:
## 1) Define variables to be read fron netCDF files
function sortVariables(defvars; onlyvars=[], addvars=[])

    outvars = Dict(:time=>"time")
    if ~isempty(onlyvars)
        # finding whether the variables are already as default:
        idx = map(x->in(x, values(defvars)), onlyvars) |> findall
        @assert ~isempty(idx) "selected only variables not present in pre-set"
        # converting only variables to Dict:
        tmp = map(x->Symbol(uppercase(x))=>x, onlyvars[idx]) |> Dict
        outvars = merge(outvars, tmp)
    else
        outvars = defvars
    end
    
    if ~isempty(addvars)
        idx = map(x->findfirst(==(x), outvars), addvars) .|> isnothing |> findall
        tmp = map(x->Symbol(uppercase(x))=>x, addvars[idx]) |> Dict
        outvars = merge(outvars, tmp)
    end
    
    return outvars
end

## 2) Read variables from netCDF file:
function retrieveVariables(ncfile::String, ncvars; attrvars=[])

    output = Dict()
    @assert isfile(ncfile) "reading $ncfile but not found!"
    ncin = NCDataset(ncfile)
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
    
    # for global attributes:
    for var ∈ attrvars
        str_var = var[2]
        println(str_var)
        tmp_var = ncin.attrib[str_var]
        # filling the output variable:
        key_var = var[1]
        output[key_var] = tmp_var
    end
    close(ncin)
    return output
end

# ********************************************
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
* :WSPD
* :WDIR
* :θ

"""
function getSondeData(sonde_file::String; addvars=[], onlyvars=[] )

    # defaul netCDF variables to read from Radiosonde:
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
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(sonde_file, ncvars)

    return output
end


# ***********************************************
# Reading Ceilomater data
"""
Function getCeil10mData(file_name::String)

This will read the ARM datafile 'file_name.nc' and return a Dictionary
with the defaul data fields.
The default data fields are:
* :time
* :heitht
* :β
* :β_raw
* :CBH
* :CBD
"""
function getCeil10mData(sonde_file::String; addvars=[], onlyvars=[], attrvars=[])

    # defaul netCDF variables to read from Ceilometer:
    ncvars = Dict(:time=>"time",
                  :height=>"range",
                  :β_raw=>"backscatter",
                  :CBH=>"first_cbh",
                  :ALT=>"alt",
                  :TILT=>"tilt_angle")

    attrib = Dict(:location=>"location_description",
                  :instrumentmodel=>"ceilometer_model",
                  :doi => "doi")
    
    # for CEIL10m data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(sonde_file, ncvars, attrvars=attrib)

    return output
end

## *******************************************************************
## Function to estimate β from raw lidar backscattering
function calculate_β_from_raw(lidar::Dict; noise_params = (100, 1e-12, 2e-7, (1.1e-8, 2.9e-8)))::Matrix{Float64}

    range_square = (lidar[:height]*1e-3).^2;  # [km²]
    beta_new = Array(1e-7*lidar[:β_raw]./range_square);  # converting to [sr⁻¹ m⁻¹] 
    signal_var = Statistics.var(beta_new[end-noise_params[1]:end,:], dims=1);
    is_saturation = findall(signal_var[:] .< noise_params[2]);
    noise = Statistics.std(beta_new[end-noise_params[1]:end,:], dims=1);
    noise_min = noise_params[3];
    noise[noise .< noise_min] .= noise_min;
    # Reset low values above Saturation:
    saturation_noise = noise_params[4][1]
    for sat_prof ∈ is_saturation
        profile = beta_new[:, sat_prof]
        peak_ind = argmax(profile)
        alt_ind = findall(profile[peak_ind:end] .< saturation_noise) .+ peak_ind .- 1
        beta_new[alt_ind, sat_prof] .= NaN
    end
    snr = beta_new./noise;
    snr_limit = 5
    β = Array(beta_new)
    ind_snr_limit = findall(snr .< snr_limit)
    β[ind_snr_limit] .= NaN
    
    return β
end

end # module
# Main file containing the package module
# See LICENSE
