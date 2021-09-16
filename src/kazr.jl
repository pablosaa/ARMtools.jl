# Part of ARMtools.jl
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
* 
* 
* lat => North latitude
* lon => East longitude
* alt => altitude above mean sea level

Attributes:
* radar_operating_frequency_chrip
* location_description
* process_version
"""
function getKAZRData(input_file::String; addvars=[], onlyvars=[], attrvars=[])
    # defaul netCDF variables to read from KAZR:
    ncvars = Dict(:time=>"time",
                  :lat=>"lat",
                  :lon=>"lon",
                  :alt=>"alt")

    attrib = Dict()
    if isvariablein(input_file, "height")
        # ARSCL file
        merge!(ncvars, Dict(:height=>"height",
                            :Ze=>"reflectivity",
                            :MDV=>"mean_doppler_velocity",
                            :SPW=>"spectral_width",
                            :LDR=>"linear_depolarization_ratio",
                            :SNR=>"signal_to_noise_ratio",
                            )
               )
        
        merge!(attrib, Dict(:location=>"location_description",
                            :instrumentmodel=>"process_version",
                            :radar_frequency => "radar_operating_frequency_chirp",
                            :doi => "doi"
                            )
               )
    elseif isvariablein(input_file, "reflectivity_copol")
        # for KAZR GE or MD data file:
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
    else
        @error "Radar data $input_file not supported...!"
    end
    
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = retrieveVariables(input_file, ncvars, attrvars=attrib)

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
        velocity_bins = range(-velocity_bins, velocity_bins, length=spectrum_n_samples)
        output[:vel_nn] = collect(velocity_bins)
    end

    # converting Attribute variables into numeric values:
    
    
    return output;
end
# ----/


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
		S = 10f0.^(0.1Zη[:, idx])
		
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
		Inoise2, Pnoise2 = findfirst(R₂ .≤ 1) |> x->!isnothing(x) ? (x, P[x]) : (10, Sₙ[10])
		
		Qnoise2 = Q[Inoise2]
				
		# For R₁ :
		Inoise1, Pnoise1 = findlast(0 .≤ diff(R₁) .< .01) |> x->!isnothing(x) ? (x, P[x]) : (10,Sₙ[10])
		Qnoise1 = Q[Inoise1]
				
		# Correcting spectral reflectivity from noise level:
		Zη_cor[:, idx], SNR_dB[:, idx] = let tmp = (S .- Pnoise2)
			tmp[tmp.≤0] .= NaN
			ii = tmp .≤ Pnoise2 + sqrt(Qnoise2)
			NoiseMean = mean(tmp[ii])	
			NoisePeak = maximum(tmp[ii])
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

# ********************************************************
# Integrate spectral reflectivity
function ∫zdη(η::T; i₀=1, i₁=size(η, 1)) where T<:AbstractArray
    i₀ = max(i₀, 1)
    i₁ = min(i₁, length(η))
    Znn = let ζnnn = @. 10f0^(0.1η)
        sum( ζnn[i₀:i₁] , dims=1)
    end
    return 10log10.(Znn)
end
# ----/

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
