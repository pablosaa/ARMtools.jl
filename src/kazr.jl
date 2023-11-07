# Part of ARMtools.jl
# sub file including functions and routines for the KAZR cloud radar

# ***************************************************************************
# Function to read KAZR ARSCRL data
"""
Function getKAZRData(file\\_name::String)

```julia-repl
julia> radar = getKAZRData("kazr_datafile.nc")
julia> radar = getKAZRData("kazr_datafile.nc", snr_filter=25)
julia> radar = getKAZRData("kazr_datafile.nc", snr_filter=nothing)
```
This will read the ARM datafile 'kazr\\_datafile.nc' and return a Dictionary
with the defaul data fields.
The default data fields are:
* :time
* :heitht
* :Ze
* :MDV
* :SPW
* :LDR
* :SNR
* :RR

OPTIONAL:
* snr_filter (Default 20%): to filter data using the percentage of (SNR\\_max - SNR\\_min)

Adiotional variables in the NetCDF file can be read by:
* addvars = ["lat", "lon", "alt"] , default []

Selected varaibles can be read by:
* onlyvars = ["Ze"], default []

Selected attributes can be read by:
* attrvars = ["process\\_version"], default []

Example of Attributes:
* radar\\_operating\\_frequency\\_chrip
* location\\_description
* process\\_version
"""
function getKAZRData(input_file::String; addvars=[], onlyvars=[], attrvars=[],
                     snr_filter=20, extras...)
    
    # defaul netCDF variables to read from KAZR:
    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    attrib = Dict()

    # checking which data level the input file is:
    data_level = giveme_datalevel(input_file)
    
    if isvariablein(input_file, "height")
        # ARSCL file
        merge!(ncvars, Dict(:height=>"height",
                            :Ze=>"reflectivity",
                            :MDV=>"mean_doppler_velocity",
                            :SPW=>"spectral_width",
                            :LDR=>"linear_depolarization_ratio",
                            :SNR=>"signal_to_noise_ratio",
                            :RR=>"precip_mean",
                            )
               )
        
        merge!(attrib, Dict(:location=>"location_description",
                            :instrumentmodel=>"process_version",
                            :radar_frequency => data_level=="c0" ? "radar_operating_frequency_chirp" : "radar_operating_frequency",
                            :doi => "doi"
                            )
               )
    elseif isvariablein(input_file, "reflectivity_copol")
        # for KAZR GE or MD data file:
        #= NOTE: KAZR GE variable names and units, e.g. for frequency depends on type of product datastream:
        for instance mosaic datastream="moskazrcfrgeM1.a1" frequency is variable and in Hz
        whereas for nsa datastream="nsakazrgeC1.a1" frequency is a global attribute :radar_operating_frequency = "34.830000 GHz" ;
        So, the assignation of variables to read needs to be done not onyl based on product level "a1", "c1", etc. but also
        based on the datastream like "KAZR GE" or "KAZR CFR GE".
        For ARSCL level c0, the frequency is a global attribute :radar_operating_frequency_chirp = "34.890000 GHz"
        This correction need to be implemented. 
        =#
        merge!(ncvars, Dict(:height => "range",
                            :Ze => "reflectivity_copol",
                            :MDV => "mean_doppler_velocity_copol",
                            :SPW => "spectral_width_copol",
                            :Zxpol => "reflectivity_xpol",  # temporal
                            :SNR => "signal_to_noise_ratio_copol",
                            )
               )
        merge!(attrib, Dict(:location=>"facility_id",
                            :instrumentmodel=>"process_version",
                            :radar_frequency => "radar_operating_frequency",
                            :fft_len => "fft_len",
                            :nyquist_velocity => "nyquist_velocity",
                            :number_spectral_ave => "num_spectral_averages",
                            :prf => "pulse_repetition_frequency",
                            :drg => "range_gate_spacing",
                            )
               )
    elseif data_level=="a1" #isinstrument(input_file, "MWACR")
        # for MWACR type radar:
        merge!(ncvars, Dict(:lat => "latitude",
                            :lon => "longitude",
                            :alt => "altitude_agl",
                            :height => "range",
                            :Ze => "reflectivity",
                            :MDV => "mean_doppler_velocity",
                            :SPW => "spectral_width",
                            :Zxpol => "reflectivity_crosspolar_v",  # temporal
                            :SNR => "signal_to_noise_ratio_copolar_h",
                            :radar_frequency => "frequency",
                            )
               )
        merge!(attrib, Dict(:location=>"facility_id",
                            :instrumentmodel=>"instrument_name",
                            :fft_len => "fft_len",
                            :nyquist_velocity => "nyquist_velocity",
                            :number_spectral_ave => "num_spectral_averages",
                            :prf => "pulse_repetition_frequency",
                            :drg => "range_gate_spacing_m",
                            )
               )
    else
        @error "Radar data $input_file not supported...!"
    end
    
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(input_file, ncvars, attrvars=attrib)

    # =========================================================
    # filtering low reflectivity values if snr_filter is TRUE:
    if !isnothing(snr_filter) && haskey(output, :SNR)
        !(0 < snr_filter < 100) && @error("snr_filter optional value out of range: ")
        filter_by_snr_threshold(output, snr_lim=snr_filter) 
    end
    
    return output
end
# ----/
function getKAZRData(rad_file::Vector{String}; addvars=[], onlyvars=[], attrvars=[], snr_filter=20, extras...)
    rad_out = Dict{Symbol, Any}()
    catvar = Dict{Symbol, Union{Nothing, Int}}()
    ntime = -1

    getdim(x,n) = findall(==(n), size(x))
    
    foreach(rad_file) do fn
        # reading single file:
        rad = getKAZRData(fn,
                          addvars=addvars,
                          onlyvars=onlyvars,
                          attrvars=attrvars,
                          snr_filter=snr_filter, extras...)

        if isempty(rad_out)
            rad_out = rad
            catvar = let ntime = length(rad[:time])
                tmp = [k=>getdim(v, ntime) for (k,v) ∈ rad if typeof(v)<:Array]
                filter(p->!isnothing(p.second), tmp) |> Dict
            end
        else
            # merging every variable following time dimension v
            [rad_out[k] = cat(rad_out[k], rad[k]; dims=v) for (k,v) ∈ catvar if !isempty(v)]
        end
    end
    return rad_out
end
# ----/


# *************************************************************************************
# Function to estimate the Doppler velocity vector from given Nyquist velocty
# and number of spectral bins.
#
"""
Return the Doppler Velocity given the Nyquist velocity and the number of spectral bins
Vnn = DopplerVelocityVector(Vn, Nspc)

Vnn is a vector of length Nspc and the zero velocity is located at index: Nspc/2 +1

"""
function DopplerVelocityVector(Vn::Real, Nspc::Int)
    # index for the location of zero Doppler velocity:
    N₀ = Nspc/2
    
    # velocity resolution of Doppler spectrum:
    Δv = Vn/N₀

    # returning vector of Doppler velocity:
    return range(-Vn, stop=Vn-Δv, length=Nspc) |> collect
end
# ----/


# ***************************************************************************
# Function to read KAZR SPECCOPOL data
"""
Function readSPECCOPOL(radar_file::String)
Function readSPECCOPOL(radar_file::String; addvars=[], onlyvars=[], attrvars=[])
Function readSPECCOPOL(radar_file::Vector{String}; addvars=[], onlyvars=[], attrvars=[])

Routines to process KAZR copol Doppler spectrum.
The default data fields are:
* :time   | Vector{DateTime}
* :height | Vector{Float} [m]
* :η_hh   | Matrix{Float} [dBm]
* :nyquist_vel | Float [m/s]
* :vel_nn | Vector{Float} [m/s]

Alternative variables can be:
* lat => North latitude
* lon => East longitude
* alt => altitude above mean sea level

Attributes:
* side_id
* doi

Note: in case radar_file is a Vector{String}, the user needs to ensure that the
file names are sorted in such a way the time dimension monotonically increases,
otherwise the time dimension will be mixed up.

"""
function readSPECCOPOLc0(kazr_file::String; addvars=[], onlyvars=[], attvars=[], extras...)

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
function readSPECCOPOL(kazr_file::Vector{String}; addvars=[], onlyvars=[], attvars=[], extras...)
    # for a list of files:
    spec_out = Dict{Symbol, Any}()
    catvar = Dict(:time=>1, :η_hh=>2, :spect_mask=>2)

    foreach(kazr_file) do fn
        spec = readSPECCOPOL(fn,
                             addvars=addvars,
                             onlyvars=onlyvars,
                             attvars=attvars, extras...)

        if isempty(spec_out) #fn==kazr_file[1]
            spec_out = spec
        else
            spec[:spect_mask][spec[:spect_mask] .≥ 0] .+= size(spec_out[:η_hh], 2)
            foreach(catvar) do (k,v)
                spec_out[k] = cat(spec_out[k], spec[k]; dims=v)
            end
        end
    end
    return spec_out
end

function readSPECCOPOL(kazr_file::String; addvars=[], onlyvars=[], attvars=[], extras...)

    # Checking which spectrum data have the file:
    if isvariablein(kazr_file, "spectra")
        ncvars = Dict(:time=>"time_offset",
                      :height=>"range",  # [m]
                      :η_hh=>"spectra", #[dbm]
                      :vel_nn=>"velocity_bins",  # [m/s]
                      :spect_mask=>"locator_mask")

        attvars = Dict(:nyquist_vel => "nyquist_velocity",
                       :calib_const => "cal_constant",
                       :Nfft => "fft_len",
                       :noise => "rx_noise",
                       :frequency => "radar_operating_frequency",
                       :Nspec_ave => "num_spectral_averages",
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
        # velocity_bins = range(-velocity_bins, velocity_bins, length=spectrum_n_samples)
        output[:vel_nn] = DopplerVelocityVector(velocity_bins, spectrum_n_samples)  #collect(velocity_bins)
    end

    # converting Attribute variables into numeric values:
    
    
    return output;
end
# ----/


# ********************************************************
# RADAR HELPER FUNCTIONS:
# ********************************************************
"""
Function return the bool Matrix with true elements of matrix fulfills
the condition that they are below the p % of the range [min max] of Matrix.
```julia-repl
julia> SNR_below_20percent = flag_array_below_lim(SNR, p=20)
```
"""
function flag_array_below_lim(X; p=20)
    # converting % to real number:
    thr0=1f-2p
    # returning bool matrix:
    X .< (filter(!isnan, X) |> extrema |> z->thr0*(z[2]-z[1]) + z[1])
end
# ----/

"""
Function to filter radar data based on a threshold of SNR:
```julia-repl
julia> filter_by_snr_threshold(radar)
```
Note that radar::Dict should already include the symbol :SNR containing the
matrix with signal-to-noise-ratio to use for the filtering.

To filter only one radar variable e.g. Ze use:
```julia-repl
julia> filter_by_snr_threshold(radar, vars=(:Ze,))
```
"""
function filter_by_snr_threshold(data::Dict; snr_lim=20, vars=())
        
    bot_val = nothing
    snr_dims = ()

    if haskey(data, :SNR)
        snr_dims = size(data[:SNR])
        bot_val = flag_array_below_lim(data[:SNR], p=snr_lim)
    else
        @warn "input data does not contain :SNR "
        return nothing
    end
    
    for (k,V) ∈ data
        k == :SNR && continue
        !isempty(vars) && all(k .!= vars) && continue
        !(typeof(V) <: Array) && continue
        (size(V) != snr_dims) && continue
        data[k][bot_val] .= NaN
    end

    return true
end
# ----/

# ********************************************************
# SPECTRUM HELPER FUNCTIONS:
# ********************************************************
# *************************************************************************************
# Following Hildebrand P.H. and Sekhon R.S. (1974) Paper:
# "Objective Determination of the Noise Level in Doppler Spectra"
# Journal of Applied Meteorology, Vol 13.
#
function Extract_Spectra_NL(Zη::Matrix; p::Int=1)
    #p = Zη[:Nspec_ave][1]
	
    # Defining output variables:
    Zη_cor = similar(Zη)
    SNR_dB = similar(Zη)
	
    Nspec, Nsamples = size(Zη)
    for idx ∈ 1:Nsamples
	# spectral reflectivity in linear units:
	S = 10f0.^(0.1Zη[:, idx]) .|> Float32
		
	# sorting S from smallest to largest:
	Iₙ = sortperm(S)
	Sₙ = S[Iₙ]
		
	ΣSₙ = cumsum(Sₙ, dims=1)
	N = cumsum(isfinite.(Sₙ), dims=1)
	fₙ = Iₙ./Nspec
		
	# eq. (4) :
	σ² = cumsum(Sₙ.*fₙ.^2, dims=1)./ΣSₙ .- (cumsum(Sₙ.*fₙ,dims=1)./ΣSₙ).^2
		
	# eq. (5)
	σₙ² = 1/12 #(F.^2)/12
		
	# eq. (6) :
	P = ΣSₙ./N
		
	# eq. (7) :
	Q = cumsum(Sₙ.^2, dims=1)./N .- P.^2
		
	# eq. (8) :
	R₁ = σₙ²./σ²
		
	# eq. (9) :
	R₂ = (P.^2)./(Q*p)
		
	# Applying criteria to determine Noise Levels
	# * For R₂ :
	Inoise2, Pnoise2 = findfirst(R₂ .≤ 1) |> x->!isnothing(x) ? (x, P[x]) : (5, Sₙ[5])
	
	Qnoise2 = Q[Inoise2] |> abs
				
	# For R₁ :
	Inoise1, Pnoise1 = findlast(0 .≤ diff(R₁) .< .01) |> x->!isnothing(x) ? (x, P[x]) : (5, Sₙ[5])
	Qnoise1 = Q[Inoise1]
				
	# Correcting spectral reflectivity from noise level:
	Zη_cor[:, idx], SNR_dB[:, idx] = let tmp = (S .- Pnoise2)
	    tmp[tmp.≤0] .= NaN
            ii = findall(tmp .≤ (Pnoise2 + sqrt(Qnoise2))) #tmp .≤ (Pnoise2 + sqrt(Qnoise2))
	    NoiseMean = mean(tmp[ii])	
            NoisePeak = isempty(ii) ? Pnoise2 : maximum(tmp[ii])
	    NoiseStdv = std(tmp[ii])
	    tmp[(tmp .≤ NoisePeak) .| isnan.(tmp)] .= Pnoise2
			
	    # Signal-to-noise-ratio:
	    SNR = tmp./Pnoise2
			
	    10log10.(tmp), 10log10.(SNR)
	end
    end
    return Zη_cor, SNR_dB
end
# ----/

# Integrate spectral reflectivity
function ∫zdη(η::T; i₀=1, i₁=size(η, 1)) where T<:AbstractArray
    i₀ = max(i₀, 1)
    i₁ = min(i₁, length(η))
    Znn = let ζnn = @. 10f0^(0.1η)
        sum( ζnn[i₀:i₁] , dims=1)
    end
    return 10f0log10.(Znn)
end
# ----/

# ********************************************************
# extract 2D Spectrogram given a time index:
function extract2DSpectrogram(spec::Dict, idx::Int; var::Symbol=:η_hh, fill_val=NaN)
    # obtaining number of samples and altitudes:
    Nspc = length(spec[:vel_nn])
    Nrng = length(spec[:height])
    # extracting indexes for spectogram height vs doppler velocity:
    idx_alt = spec[:spect_mask][:, idx] .≥ 0
    idx_rng = spec[:spect_mask][idx_alt, idx] .+ 1
    # creating output 2D variable:
    out2D = fill(fill_val, Nrng, Nspc)
    out2D[idx_alt, :] = spec[var][:, idx_rng]'

    rng_lim = findall(!isnan, out2D[:,1]) |> extrema |> x-> spec[:height][[x...]] .|> floor
	
    return out2D, rng_lim
end
# ----/

#= ********************************************************
Function to reshape matrix or vector M into a 2D shape given by "masked"
masked needs to be a bool matrix or Int matrix
=#
function maskedReshape(M::T, masked::T; Dim3=nothing) where T<:AbstractArray
    # verifying dimensions match:
    (prod∘sum)(M) == Dim3*length(masked[:])
    xy = size(masked)
    DimXY = isnothing(Dim3) ? xy : (xy...,Dim3)
    M_out = fill(NaN32, DimXY)
    idx_xy = cumsum(masked[:]) |> A->reshape(A, xy)
    for i ∈ 1:xy[1]
        for j ∈ 1:xy[2]
            !masked[i,j] && continue
            if isnothing(Dim3)
                M_out[i,j] = M[idx_xy[i,j]]
            else
                M_out[i,j,:] = M[idx_xy[i,j],:]
            end
        end
    end
    return M_out
end
# ---/

# end of script.
# *****************************************************************
##function η_dBz(spec::Dict)
##
##    # calculating spectral reflectivity:
##    (nh, nt) = size(spec[:spect_mask])
##    Pₙ = spec[:noise][1]    # dBm
##    C  = spec[:cal_constant][1]   # dB
##    η_out = similar(spec[:η_hh])
##    
##    for it ∈ (1:nt)
##        for ih ∈ (1:nh)
##            k = spec[:spect_mask][ih, it]
##            if k ≥ 0
##                k += 1
##            else
##                continue
##            end
##            
##            let range_factor = 20log10(1f-3spec[:height][ih]),
##                    Pₛ = spec[:η_hh][:, k]  # dBm
##            
##                    Pₜ = Pₛ .- Pₙ    # dB
##                    Znn = @. Pₜ + range_factor - C # dB
##                    η_out[:,k] = denoise(Znn, wavelet(WT.sym6, WT.Filter), dnt=VisuShrink(256), TI=true)
##                
##            end
##        end
##    end
##    return η_out
##end
##function η_dBz(filen::String)
##    @assert isfile(filen) "$filen seems not to exist!"
##    spec = ARMtools.readSPECCOPOL(filen)
##    spec[:Zη] = η_dBz(spec)
##    return spec
##end
### ----/
