using Images, ColorSchemes

function imwrite(filename, M; cmap=ColorSchemes.seismic, lims=0)

    N = length(cmap)

    # Get array limits:

    if lims==0
        lims = [minimum(M), maximum(M)]
    end

    M[M .< lims[1]] .= lims[1]
    M[M .> lims[2]] .= lims[2]

    # Rescale array:

    M = (N - 1) / (lims[2] - lims[1]) * (M .- lims[1]) .+ 1
    M = Int.(round.(M)) 

    # Create image:

    I = [cmap[i] for i in M']

    # Save image:

    save(filename, I)

    return nothing

end

function save_frames(framename, M; lims=0)

    n = size(M)[3]

    if lims == 0
        l = max(maximum(M), -minimum(M))
        lims = [-l, l]
    end

    for i in 1:n

        framenum = lpad(string(i), Int(floor(log10(n))+1), '0')
        filename = framename * "_" * framenum * ".png"

        imwrite(filename, M[:, :, i]; cmap=ColorSchemes.seismic, lims=lims)

    end

    return nothing

end

nothing