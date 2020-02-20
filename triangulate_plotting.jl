# ## 0 Preliminaries

import Triangulate
# using PyPlot
using Plots
# using ExtendableSparse
# using SparseArrays
using Printf
# using BenchmarkTools
# using DataFrames

# These backends have been tested:
gr() # 👍
#pyplot() # 👍 - but built-in method should be used
#unicodeplots() # 👎 - UnicodePlots does not handle line segments correctly
#plotly() # 👍 - but can be improved

plot(size=(500,400), aspect_ratio=:equal, legend=false)

# # From SciComp homework:

"Create coordinates of `k` vertices on the circle in polar coordinates"
function circleapprox(k::Int)
    phis = LinRange(0, 2 * pi, k + 1)[1:end-1]

    x = map(cos, phis)
    y = map(sin, phis)

    pointlist = [x y]

    segmentlist = [collect(1:k) collect(2:k+1)]
    segmentlist[end, end] = 1

    segmentmarkerlist = ones(Int64, k)

    return [pointlist, segmentlist, segmentmarkerlist]
end


function task1(
    ;
    plotgrid = false,
    nref0 = 1,
    nref1 = 7,
    k = 1,
    l = 1,
    N = 36,
    verbose = false,
)
    polygon = circleapprox(N)
    if nref0 > nref1
        nref1 = nref0
    end
    trioutlist = []


    for iref = nref0:nref1
        triin = Triangulate.TriangulateIO()
        triin.pointlist = Matrix{Cdouble}(polygon[1]')
        triin.segmentlist = Matrix{Cint}(polygon[2]')
        triin.segmentmarkerlist = Vector{Int32}(polygon[3])
        area = 0.1 * 2.0^(-2 * iref)
        s_area = @sprintf("%.20f", area)


        (triout, vorout) = Triangulate.triangulate("pqa$(s_area)qQD", triin)
        n = size(triout.pointlist, 2)

        trioutlist = push!(trioutlist, (triout, vorout))
        ntest = size(trioutlist[iref][1].pointlist, 2)
    end

    return trioutlist
end



meshes = task1(plotgrid = false, verbose = false, nref0=1, nref1=1)
triout, vorout = meshes[1]

exact_solution(x, y) = sin(x*pi) + sin(y*pi)

# Discretize
u_discretised = map(1:size(triout.pointlist, 2)) do i
    return exact_solution(
        triout.pointlist[1, i],
        triout.pointlist[2, i],
    )
end

function tripcolor(triout, point_values)

    # We add a NaN point
    tp = triout.pointlist
    p = [triout.pointlist [NaN; NaN]]

    nanindex = size(p, 2)


    trilist = triout.trianglelist[:,:]
    trilist_extended = [trilist; trilist[[1],:]; fill(nanindex, (1, size(trilist, 2)))]

    # Draw triangle mesh
    # plot!(p[1, trilist_extended[:]], p[2, trilist_extended[:]], linecolor="black", linewidth=1, linealpha=.2)

    # Scatter values at points
    # plot!(tp[1,:], tp[2,:], marker_z=point_values, markerstrokewidth=0, markersize=5,seriestype=:scatter)

    tri_colors = Array{Float64, 1}(undef, size(trilist, 2))

    for tri_index in 1:size(trilist, 2)
        p_indices = trilist[:,tri_index]
        points = tp[:,p_indices]

        values = point_values[p_indices]

        tri_colors[tri_index] = sum(values) / 3.0
    end

    plot!(p[1, trilist_extended[:]], p[2, trilist_extended[:]], linecolor="black", linewidth=1, linealpha=.2, seriestype=:shape, fill_z=tri_colors, fillcolor=:grays)
end

tripcolor(triout, u_discretised)

##

function low_mid_high(contestors, vals)
    # Equivalent to: sort(vals[contestors] |> enumerate |> collect, by=t->t[2])
    a, b, c = vals[contestors]

    if a <= b
        if b <= c
            return (1, a), (2, b), (3, c)
        else
            if a <= c
                return (1, a), (3, c), (2, b)
            else
                return (3, c), (1, a), (2, b)
            end
        end
    else
        # b < a
        if a <= c
            return (2, b), (1, a), (3, c)
        else
            # b < a; c < a
            if b <= c
                return (2, b), (3, c), (1, a)
            else
                return (3, c), (2, b), (1, a)
            end
        end
    end
end


function interpolate(two_points :: Array{Float64, 2}, t_bounds, t_value)
    # `two_points` must have the points as _columns_
    α = (t_value - t_bounds[1]) / (t_bounds[2] - t_bounds[1])
    return (1.0 - α) * two_points[:, 1] .+ α * two_points[:, 2]
end


function tricontour(triout, point_values; num_contours = 11)
    val_min, val_max = extrema(point_values)

    # We transform all point values into the range [0, num_contours], and then round down.
    # These integer values give the index of the contour outside of that point.
    u_rounded = Array{Int64}(floor.(num_contours .* (point_values .- val_min) ./ (val_max - val_min) ))

    # (takes about 10% of runtime)


    contour_paths = Array{Array{Float64,1},1}()
    contour_path_colors = Array{Float64,1}()

    for tri_index in 1:size(triout.trianglelist, 2)
        p_indices = view(triout.trianglelist,:,tri_index)
        points = view(triout.pointlist, :,p_indices)

        values = point_values[p_indices]
        values_rounded = u_rounded[p_indices]

        # low, mid and high are the indices of the three triangle verts, sorted by rounded value.
        # (low, low_v), (mid, mid_v), (high, high_v) = sort(values_rounded |> enumerate |> collect, by=t->t[2])
        (low, low_v), (mid, mid_v), (high, high_v) = low_mid_high(p_indices, u_rounded)

        # We focus on the edge between the `low` vertex and the `high` vertex, and iterate over intermediate contour values. (All contours that intersect the triangle must intersect this edge.)
        # (If no contour lines pass through the triangle, then `low_v == mid_v == high_v`, and the iteration is empty.)
        for contour_value in (low_v + 1):high_v
            denormalised_value = val_min + contour_value * (val_max - val_min) / num_contours

            # We push the color value twice, because the line segment has a start coordinate, end coordinate and NaN coordinate.
            push!(contour_path_colors, denormalised_value, denormalised_value, denormalised_value)

            # On the `low` - `high` edge, we find the contour intersection by linear interpolation:
            push!(contour_paths, interpolate(points[:,[low, high]], values[[low, high]], denormalised_value))

            # Will this contour intersect the `low` - `mid` edge or the `mid` - `high` edge?
            if contour_value <= mid_v
                push!(contour_paths, interpolate(points[:,[low, mid]], values[[low, mid]], denormalised_value))
            else
                push!(contour_paths, interpolate(points[:,[mid, high]], values[[mid, high]], denormalised_value))
            end

            # push a NaN coordinate to terminate the path.
            push!(contour_paths, [NaN, NaN])
        end
    end

    # push!(transitions, [1.0, 1.0])

    # Reformat into plottable structure:
    path_array = hcat(contour_paths...)
    plot!(path_array[1,:], path_array[2,:], line_z = contour_path_colors, linewidth=5)
end

tricontour(triout, u_discretised)

savefig("triangulate_plotting_demo_" + backend_name() + ".pdf")
savefig("triangulate_plotting_demo_" + backend_name() + ".png")

plot!()
