# Part of ARNtools.jl
# containing functions related to NAV products e.g. RV Polarstern
# during MOSAiC

# TODO: optional time steps to read line Minute, 30Seconds, etc.
# *************************************************
# Function to read NAV products
function getNAVData(in_file::String; addvars=[], onlyvars=[], attrvars=[],
                    time_res=Dates.Minute(1))
    ncvars = Dict(:time => "time",
                  :lat => "lat",
                  :lon => "lon",
                  :alt => "alt",
                  :pitch => "pitch",
                  :heave => "heave",
                  :yaw => "yaw",
                  :roll => "roll")

    # for NAV data
    ncvars = sortVariables(ncvars, onlyvars=onlyvars, addvars=addvars)

    output = Dict()
    let dummy = retrieveVariables(in_file, ncvars, attrvars=attrvars)
        idx_ts = index_at_time_resolution(dummy[:time], δts=time_res)
        foreach(keys(ncvars) |> collect) do kk
            output[kk] = dummy[kk][idx_ts]
        end
    end
    
    return output
end


# end of script
