"""
Julia Package including functions to read and process data from the DOE ARM program.

The package is mainly focused to work with instruments from the NSA facility in Utqiagvik, Alaska.
Nonetheless since most of the ARM instrumentation is standard, the package can also be used for
other facilities.

(c) 2020, Pablo Saavedra Garfias
University of Leipzig
Faculty of Physics and Geosciences
LIM

See LICENSE
"""
module ARMtools

using NCDatasets, Statistics, Printf

## * Auxiliary functions:
## 1) Define variables to be read fron netCDF files
function sortVariables(defvars; onlyvars=[], addvars=[])

    outvars = Dict(:time=>"time")
    if ~isempty(onlyvars)
        # finding whether the variables are already as default:
        idx = map(x->in(x, onlyvars), values(defvars)) |> findall
        @assert ~isempty(idx) "selected only variables not present in pre-set"

        # converting only variables to Dict:
        thekeys = collect(keys(defvars))
        thevalues = collect(values(defvars))
        tmp = map(x->thekeys[x] => thevalues[x], idx) |> Dict
        #tmp = map(x->Symbol(uppercase(x))=>x, onlyvars[idx]) |> Dict
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
# ----/

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
# ----/

# ****************************************
# * get file from pattern:
"""
 Function getFilePattern(path::String, product::String, yy, mm ,dd)

 retrieve file name pattern based on year, month, day to get
 a string to read.
"""
function getFilePattern(path::String, product::String, yy, mm ,dd)
    base_dir = joinpath(path, product, @sprintf("%04d", yy))
    list_file = readdir(base_dir, join=true)
    pattern = @sprintf("%04d%02d%02d", yy, mm, dd)
    ofile = filter(x->all(occursin.(pattern, x)), list_file)
    ofile = isempty(ofile) ? "$dd.$mm.$yy none" : ofile[1]  
    return ofile
end
# ----/


# ********************************************
# ********************************************
# READING FUNCTIONS

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
# ----/

# ***********************************************
# Reading data for Ceilomater type Vaisala 30
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
# ----/

# ***********************************************
# Reading data for HSRL lidar
"""
Function getHSRLData(file_name::String)

This will read the ARM datafile 'file_name.nc' and return a Dictionary
with the defaul data fields.
The default data fields are:
* :time
* :heitht
* :β
* :β_raw
* :δ
* :snr

Alternative variables can be:
- beta_a => Particulate extinction cross-section per unit volume
- atten_beta_r_backscatter => Attenuated molecular return
- lat => North latitude
- lon => East longitude
- alt => altitude above mean sea level
"""
function getHSRLData(sonde_file::String; addvars=[], onlyvars=[], attrvars=[])

    # defaul netCDF variables to read from HSRL:
    ncvars = Dict(:time=>"time",
                  :height=>"range",
                  :β_raw=>"beta_a_backscatter",
                  :SNR=>"beta_a_backscatter_std",
                  :δ=>"depol")

    attrib = Dict(:location=>"dod_version",
                  :instrumentmodel=>"facility_id")
    
    # for CEIL10m data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(sonde_file, ncvars, attrvars=attrib)

    return output
end
# ----/

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
# ----/


# ***************************************************************************
# Function to read KAZR ARSCRL data
"""
Function getKAZRData(file_name::String)

This will read the ARM datafile 'file_name.nc' and return a Dictionary
with the defaul data fields.
The default data fields are:
* :time
* :heitht
* :Ze
* :MDV
* :SPW
* :LDR
* :SNR

Alternative variables can be:
* beta_a => Particulate extinction cross-section per unit volume
* atten_beta_r_backscatter => Attenuated molecular return
* lat => North latitude
* lon => East longitude
* alt => altitude above mean sea level

Attributes:
* radar_operating_frequency_chrip
* location_description
* process_version
"""
function getKAZRData(sonde_file::String; addvars=[], onlyvars=[], attrvars=[])
    # defaul netCDF variables to read from KAZR:
    ncvars = Dict(:time=>"time",
                  :height=>"height",
                  :Ze=>"reflectivity",
                  :MDV=>"mean_doppler_velocity",
                  :SPW=>"spectral_width",
                  :LDR=>"linear_depolarization_ratio",
                  :SNR=>"signal_to_noise_ratio",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    attrib = Dict(:location=>"location_description",
                  :instrumentmodel=>"process_version")
    
    # for KAZR data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(sonde_file, ncvars, attrvars=attrib)

    return output
end
# ----/


# *************************************************
# Function to read MWR RET products
"""
Function getMWRData(mwr_file::String; addvars=[], onlyvars=[], attrvars=[])

This function read data from the Microwave Radiometer RET retrievals

The default data fields are:
* :time
* :LWP
* :IWV

Alternative variables can be:
* lat => North latitude
* lon => East longitude
* alt => altitude above mean sea level
* surface_vapor_pres => 
* surface_pres => 
* surface_rh =>
* tbsky23 => Brightness Temperature at 23GHz [K]
* tbsky31 => Brightness Temperature at 31GHz [K]

Attributes:
* side_id
* doi
"""
function getMWRData(mwr_file::String; addvars=[], onlyvars=[], attrvars=[])

    ARM_PRODUCT = "RET"
    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    if ARM_PRODUCT=="RET"
        ncvars[:LWP] = "be_lwp";  # [g/m²]
        ncvars[:IWV] = "be_pwv";  # [cm]
        
        factor_lwp = 1f0;
        factor_iwv = 997f-2;  # [cm] -> [kg m⁻²]

        ncvars[:CBH] = "cloud_base_height"   # [km] AGL
        ncvars[:CLT] = "cloud_temp"   # [K]
        ncvars[:SFT] = "surface_temp"   # [K]
        
    elseif ARM_PRODUCT=="LOS"
        ncvars[:LWP] = "liq";  # [cm]
        ncvars[:IWV] = "vap";  # [cm]
        factor_lwp = 1f4;
        factor_iwv = 997f-2;  # [cm] -> [kg m⁻²]
    else
    end
    
    # for MWR data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(mwr_file, ncvars, attrvars=attrvars)

    # converting units for LWP and IWV
    if haskey(output, :LWP)
        output[:LWP] *= factor_lwp
    end
    if haskey(output, :IWV)
        output[:IWV] *= factor_iwv
    end
    
    return output
end
# ----/

end # module
# Main file containing the package module
# See LICENSE
