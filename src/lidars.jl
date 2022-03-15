## -----------------------------------------------------------
# Collection of functions specific for ARM Lidar systems
# This file includes reading and processing functions.
# So far included:
# * Vaisala C30
# * HSRL
# * MPL
# * Jenoptik (Lufft)
## -----------------------------------------------------------

# ***********************************************
function getLidarData(lidar_file::String)
   
    # Checking which type of LIDAR data is in:
    if isvariablein(lidar_file, "backscatter")
        # then CEIL10m
        str_func = "CEIL10m"
        ex = getCeil10mData
        
    elseif isvariablein(lidar_file, "beta_a_backscatter")
        # then HSRL
        str_func = "HSRL"
        ex = getHSRLData
    elseif isvariablein(lidar_file, "backscatter_snr")
        # then MPL
        str_func = "MPL"
        ex = getMPLData
    else
        @error "$lidar_file seems not an ARM LIDAR data file!"
    end
    
    println("loading $str_func data...")
    ex = :($ex($lidar_file))
    return eval(ex)
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
                  :height => "range",
                  :β_raw => "beta_a_backscatter",
                  :SNR => "beta_a_backscatter_std",
                  :δ => "depol")

    attrib = Dict(:location => "dod_version",
                  :instrumentmodel => "facility_id",
                  :doi => "doi")
    
    # for HSRL data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(sonde_file, ncvars, attrvars=attrib)

    return output
end
# ----/


# ***********************************************
# Reading data for MPL lidar
"""
Function getMPLData(file_name::String)

This will read the ARM datafile 'file_name.nc' and return a Dictionary
with the defaul data fields.
The default data fields are:
* :time
* :heitht
* :β
* :β_raw
* :δ
* :snr
* :cbh
* :clh

Alternative variables can be:
- :cta => cloud_top_attenuation_flag
- atten_beta_r_backscatter => Attenuated molecular return
- lat => North latitude
- lon => East longitude
- alt => altitude above mean sea level
"""
function getMPLData(sonde_file::String; addvars=[], onlyvars=[], attrvars=[])

    # defaul netCDF variables to read from MPL:
    ncvars = Dict(:time => "time",
                  :height => "height",
                  :β_raw => "backscatter",
                  :SNR => "backscatter_snr",
                  :δ => "linear_depol_ratio",
                  :cbh => "cloud_base",)

    attrib = Dict(:location => "site_id",
                  :doi => "doi")
    
    # for MPL data:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(sonde_file, ncvars, attrvars=attrib)

    return output
end
# ----/

## ******************************************
#  PROCESSING FUNCTIONS
## ******************************************
function raw_to_β(β_in::Matrix, snr::Matrix, range::Vector; snr_limit = 5, δc::Matrix=Matrix{AbstractFloat}(undef, 0, 0))
    range_square = range.^2
    range_square[1] = range[1] ≈ 0 ? range[1] : range_square[1]
    β_new = β_in./range_square
    idx_snr = snr .< snr_limit
    β_new[idx_snr] .= NaN

    if !isempty(δc)
        δₗ = circular_to_linear_depol(δc)
        δₗ[idx_snr] .= NaN
        return β_new, δₗ
    end
    
    return β_new
end
# ----/

## Circular depolarization to Linear depolarization
function circular_to_linear_depol(δ_c::Matrix)::Matrix
    return @. δ_c/(2.0 + δ_c)
end
# ----/

## *******************************************************************
## Function to estimate β from raw lidar backscattering
# Note: the noise_params optional variable is general and it should
# be replaced by parameters specific for every lidar system.
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



# end of file
# ---/
