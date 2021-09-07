# sub file including functions and routines for the KAZR cloud radar

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

# ***************************************************************************
# Function to read KAZR SPECCOPOL data
"""
Function readSPECCOPOL(kazr_file::String)
Function readSPECCOPOL(kazr_file::String; addvars=[], onlyvars=[], attrvars=[])

Routines to process KAZR copol Doppler spectrum.
The default data fields are:
* :time
* :height    | [m]
* :η_hh  | [dbm]
* :nyquist_vel | [m/s]

Alternative variables can be:
* lat => North latitude
* lon => East longitude
* alt => altitude above mean sea level

Attributes:
* side_id
* doi

"""
function readSPECCOPOLc0(kazr_file::String; addvars=[], onlyvars=[], attvars=[])

    ncvars = Dict(:time=>"time",
                  :height=>"range",  # [m]
                  :η_hh=>"radar_power_spectrum_of_copolar_h", #[dbm]
                  :nyquist_vel=>"nyquist_velocity",  # [m/s]
                  :spect_idx=>"spectrum_index",
                  :noise_sky=>"radar_measured_sky_noise_h")

    # variables selection from input:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(kazr_file, ncvars, attrvars=attvars)

    return output
end
function readSPECCOPOL(kazr_file::String; addvars=[], onlyvars=[], attvars=[])

    # Checking which spectrum data have the file:
    if isvariablein(kazr_file, "spectra")
        ncvars = Dict(:time=>"time_offset",
                      :height=>"range",  # [m]
                      :η_hh=>"spectra", #[dbm]
                      :vel_nn=>"velocity_bins",  # [m/s]
                      :spect_mask=>"locator_mask")

        attvars = Dict(:nyquist_vel => "nyquist_velocity",
                       :cal_constant => "cal_constant",
                       :Nfft => "fft_len",
                       :noise => "rx_noise",
                       :frequency => "radar_operating_frequency",
        )
    elseif isvariablein(kazr_file, "radar_power_spectrum_of_copolar_h")
        ncvars = Dict(:time=>"time",
                      :height=>"range",  # [m]
                      :η_hh=>"radar_power_spectrum_of_copolar_h", #[dbm]
                      :vel_nn=>"nyquist_velocity",  # [m/s]
                      :spect_mask=>"spectrum_index")
    elseif isvariablein(kazr_file, "radar_power_spectrum_of_copolar_v")
        @info "sorry... copolar V not yet implemented"
    else
        @error "$kazr_file seems not to be an ARM Spectrum data file"
    end

    # variables selection from input:
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(kazr_file, ncvars, attrvars=attvars)

    # converting nyquist_velocity to velocity bins:
    if values(ncvars) .|> contains("nyquist_velocity") |> any
        velocity_bins = unique(output[:vel_nn])[1]
    
        spectrum_n_samples, n_samples = size(output[:η_hh])
        velocity_bins = range(-velocity_bins, velocity_bins, length=spectrum_n_samples)
        output[:vel_nn] = collect(velocity_bins)
    end

    # converting Attribute variables into numeric values:
    
    
    return output;
end
# ----/


# *************************************************************************************
#
using Wavelets
function η_dBz(spec::Dict)

    # calculating spectral reflectivity:
    (nh, nt) = size(spec[:spect_mask])
    Pₙ = spec[:noise][1]    # dBm
    C  = spec[:cal_constant][1]   # dB
    η_out = similar(spec[:η_hh])
    
    for it ∈ (1:nt)
        for ih ∈ (1:nh)
            k = spec[:spect_mask][ih, it]
            if k ≥ 0
                k += 1
            else
                continue
            end
            
            let range_factor = 20log10(1f-3spec[:height][ih]),
                    Pₛ = spec[:η_hh][:, k]  # dBm
            
                    Pₜ = Pₛ .- Pₙ    # dB
                    Znn = Pₜ .+ range_factor .- C # dB
                    η_out[:,k] = denoise(Znn, wavelet(WT.sym6, WT.Filter), dnt=VisuShrink(256), TI=true)
                
            end
        end
    end
    return η_out
end
function η_dBz(filen::String)
    spec = ARMtools.readSPEC(filen)
    return η_dBz(spec)
end
# ----/

# end of script.
