# Ecliptic coordinate-system illustration with the coordinate frames integrated
# directly into the main figure.
#
# One-time installation:
#   julia -e 'using Pkg; Pkg.add(["CairoMakie", "LaTeXStrings"])'
#
# Run:
#   julia ecliptic_coordinate_system_main_frames.jl
#
# Output:
#   ecliptic_coordinate_system_main_frames.pdf
#   ecliptic_coordinate_system_main_frames.svg

using CairoMakie
using LinearAlgebra
using LaTeXStrings

# =============================================================================
# Geometry utilities
# =============================================================================

deg(x) = π * x / 180
unit(v) = Float64.(v) / norm(v)

struct OrthoCamera
    right::Vector{Float64}
    up::Vector{Float64}
    view::Vector{Float64}
end

function camera_frame(view; world_up = [0.0, 0.0, 1.0])
    v = unit(view)
    up0 = unit(world_up)

    if abs(dot(v, up0)) > 0.98
        up0 = [0.0, 1.0, 0.0]
    end

    right = unit(cross(up0, v))
    up = unit(cross(v, right))
    return OrthoCamera(right, up, v)
end

project_point(p, camera::OrthoCamera) =
    Point2f(dot(p, camera.right), dot(p, camera.up))

project_dir(v, camera::OrthoCamera) =
    Vec2f(dot(v, camera.right), dot(v, camera.up))

depth(p, camera::OrthoCamera) = dot(p, camera.view)

circle_point(longitude; radius = 1.0) =
    radius .* [cos(longitude), sin(longitude), 0.0]

function ecliptic_point(longitude, latitude; radius = 1.0)
    return radius .* [
        cos(latitude) * cos(longitude),
        cos(latitude) * sin(longitude),
        sin(latitude),
    ]
end

function rotate_about_axis(p, axis, angle)
    a = unit(axis)
    c = cos(angle)
    s = sin(angle)
    return c .* p .+ s .* cross(a, p) .+ (1 - c) * dot(a, p) .* a
end

# =============================================================================
# Front/back visibility splitting
# =============================================================================

const GAP2 = Point2f(NaN, NaN)
isgap(p::Point2f) = isnan(p[1]) || isnan(p[2])

function add_gap!(points::Vector{Point2f})
    if !isempty(points) && !isgap(points[end])
        push!(points, GAP2)
    end
    return points
end

function projected_side(curve, camera::OrthoCamera; front = true, closed = false)
    points3 = closed ? vcat(curve, [curve[1]]) : curve
    points2 = Point2f[]

    wanted(d) = front ? d >= 0 : d < 0

    for i in 1:(length(points3) - 1)
        a = points3[i]
        b = points3[i + 1]
        da = depth(a, camera)
        db = depth(b, camera)
        a_wanted = wanted(da)
        b_wanted = wanted(db)

        if a_wanted && b_wanted
            if isempty(points2) || isgap(points2[end])
                push!(points2, project_point(a, camera))
            end
            push!(points2, project_point(b, camera))

        elseif a_wanted && !b_wanted
            if isempty(points2) || isgap(points2[end])
                push!(points2, project_point(a, camera))
            end
            τ = da / (da - db)
            crossing = (1 - τ) .* a .+ τ .* b
            push!(points2, project_point(crossing, camera))
            add_gap!(points2)

        elseif !a_wanted && b_wanted
            add_gap!(points2)
            τ = da / (da - db)
            crossing = (1 - τ) .* a .+ τ .* b
            push!(points2, project_point(crossing, camera))
            push!(points2, project_point(b, camera))

        else
            add_gap!(points2)
        end
    end

    while !isempty(points2) && isgap(points2[end])
        pop!(points2)
    end

    return points2
end

function draw_spherical_curve!(
        ax,
        curve,
        camera;
        closed = true,
        linewidth = 1.5,
        hidden_color = (:black, 0.45),
        visible_color = :black,
        hidden_linestyle = (:dash, :loose),
    )

    rear = projected_side(curve, camera; front = false, closed = closed)
    front = projected_side(curve, camera; front = true, closed = closed)

    !isempty(rear) && lines!(
        ax,
        rear;
        color = hidden_color,
        linewidth = linewidth,
        linestyle = hidden_linestyle,
    )

    !isempty(front) && lines!(
        ax,
        front;
        color = visible_color,
        linewidth = linewidth,
    )

    return nothing
end

function add_arrowhead!(ax, curve2; backstep = 10, linewidth = 1.4, color = :black)
    n = length(curve2)
    n < 2 && return nothing

    p1 = curve2[max(1, n - backstep)]
    p2 = curve2[end]

    arrows2d!(
        ax,
        [p1],
        [p2];
        argmode = :endpoint,
        color = color,
        shaftwidth = linewidth,
        tipwidth = 10,
        tiplength = 8,
    )
    return nothing
end

# =============================================================================
# Coordinate-frame axes drawn in the main figure
# =============================================================================

function draw_projected_axis!(
        ax,
        origin2::Point2f,
        direction3,
        label;
        camera,
        length = 0.26,
        linewidth = 1.0,
        color = (:black, 0.55),
        label_shift = Vec2f(0.0, 0.0),
        fontsize = 15,
    )

    d = project_dir(direction3, camera)
    norm(d) < 1e-8 && return nothing
    d = length * d / norm(d)
    endpoint = origin2 + d

    arrows2d!(
        ax,
        [origin2],
        [endpoint];
        argmode = :endpoint,
        color = color,
        shaftwidth = linewidth,
        tipwidth = 9,
        tiplength = 7,
    )

    text!(ax,
        endpoint[1] + label_shift[1],
        endpoint[2] + label_shift[2];
        text = label,
        align = (:center, :center),
        fontsize = fontsize,
        color = color,
    )

    return nothing
end

# =============================================================================
# Main illustration
# =============================================================================

function ecliptic_coordinate_figure(;
        λ = deg(100),
        β = deg(35),
        ε = deg(23.4393),
        zero_azimuth = deg(-160),
        camera_view = [-3.0, -4.0, 2.2],
        show_full_latitude_circle = true,
        show_coordinate_frames = true,
    )

    camera = camera_frame(camera_view)

    # Ecliptic frame
    x_ecl = unit(circle_point(zero_azimuth))
    z_ecl = [0.0, 0.0, 1.0]
    y_ecl = unit(cross(z_ecl, x_ecl))

    # Equatorial frame: same x-axis, rotated z and y axes.
    x_eq = x_ecl
    z_eq = rotate_about_axis(z_ecl, x_ecl, ε)
    y_eq = unit(cross(z_eq, x_eq))

    θ = collect(range(0, 2π; length = 721))[1:end-1]
    ecliptic = [circle_point(t) for t in θ]
    equator = [rotate_about_axis(p, x_ecl, ε) for p in ecliptic]

    star_longitude = zero_azimuth + λ
    foot = ecliptic_point(star_longitude, 0.0)
    north_ecliptic_pole = z_ecl

    latitude_circle = [
        cos(t) .* foot .+ sin(t) .* north_ecliptic_pole
        for t in θ
    ]

    longitude_arc = [
        circle_point(t)
        for t in range(zero_azimuth, star_longitude; length = 181)
    ]

    latitude_arc = [
        ecliptic_point(star_longitude, b)
        for b in range(0, β; length = 121)
    ]

    star = latitude_arc[end]

    fig = Figure(size = (860, 860), figure_padding = 10)
    ax = Axis(fig[1, 1]; aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    xlims!(ax, -1.14, 1.14)
    ylims!(ax, -1.14, 1.14)

    O = Point2f(0, 0)
    p_north = project_point(north_ecliptic_pole, camera)

    # Coordinate frames inside the main sphere.
    if show_coordinate_frames
        # Common x-axis: vernal-equinox direction.
        draw_projected_axis!(
            ax, O, x_ecl, L"x_{\mathrm{ecl}} \equiv x_{\mathrm{eq}}";
            camera = camera,
            length = 1.05,
            linewidth = 0.95,
            color = (:black, 0.58),
            label_shift = Vec2f(-0.02, -0.04),
            fontsize = 14,
        )

        draw_projected_axis!(
            ax, O, y_ecl, L"y_{\mathrm{ecl}}";
            camera = camera,
            length = 1.05,
            linewidth = 0.9,
            color = (:black, 0.45),
            label_shift = Vec2f(0.01, 0.025),
            fontsize = 14,
        )

        draw_projected_axis!(
            ax, O, z_ecl, L"z_{\mathrm{ecl}}";
            camera = camera,
            length = 1.05,
            linewidth = 0.95,
            color = (:black, 0.45),
            label_shift = Vec2f(0.03, 0.02),
            fontsize = 14,
        )

        draw_projected_axis!(
            ax, O, y_eq, L"y_{\mathrm{eq}}";
            camera = camera,
            length = 1.05,
            linewidth = 0.85,
            color = (:black, 0.30),
            label_shift = Vec2f(-0.01, 0.02),
            fontsize = 14,
        )

        draw_projected_axis!(
            ax, O, z_eq, L"z_{\mathrm{eq}}";
            camera = camera,
            length = 1,
            linewidth = 0.9,
            color = (:black, 0.30),
            label_shift = Vec2f(0.02, 0.025),
            fontsize = 14,
        )
    else
        lines!(ax, [O, p_north]; color = (:black, 0.48), linewidth = 1.0)
    end

    draw_spherical_curve!(
        ax,
        equator,
        camera;
        linewidth = 1.35,
        hidden_color = (:black, 0.34),
        visible_color = (:black, 0.82),
    )

    draw_spherical_curve!(
        ax,
        ecliptic,
        camera;
        linewidth = 1.65,
        hidden_color = (:black, 0.42),
        visible_color = :black,
    )

    if show_full_latitude_circle
        draw_spherical_curve!(
            ax,
            latitude_circle,
            camera;
            linewidth = 1.05,
            hidden_color = (:black, 0.22),
            visible_color = (:black, 0.48),
        )
    end

    longitude2 = [project_point(p, camera) for p in longitude_arc]
    latitude2 = [project_point(p, camera) for p in latitude_arc]

    lines!(ax, longitude2; color = :black, linewidth = 2.15)
    add_arrowhead!(ax, longitude2; backstep = 12, linewidth = 1.45)

    lines!(ax, latitude2; color = :black, linewidth = 2.0)
    add_arrowhead!(ax, latitude2; backstep = 12, linewidth = 1.4)

    p_star = project_point(star, camera)
    arrows2d!(
        ax,
        [O],
        [p_star];
        argmode = :endpoint,
        color = :black,
        shaftwidth = 1.35,
        tipwidth = 11,
        tiplength = 9,
    )

    outline = [Point2f(cos(t), sin(t)) for t in range(0, 2π; length = 721)]
    lines!(ax, outline; color = :black, linewidth = 2.05)

    p_zero = project_point(x_ecl, camera)
    p_foot = project_point(foot, camera)

    scatter!(ax, [O]; color = :black, markersize = 5)
    scatter!(ax, [p_zero]; color = :black, markersize = 5)
    scatter!(ax, [p_foot]; color = :black, markersize = 6)
    scatter!(ax, [p_star]; color = :black, markersize = 7)

    λ_mid = project_point(circle_point(zero_azimuth + 0.55λ), camera)
    β_mid = project_point(ecliptic_point(star_longitude, 0.55β), camera)

    ecliptic_label_point = project_point(
        circle_point(zero_azimuth + deg(20)),
        camera,
    )

    equator_label_point = project_point(
        rotate_about_axis(
            circle_point(zero_azimuth + deg(72)),
            x_ecl,
            ε,
        ),
        camera,
    )

    text!(ax, p_north[1], p_north[2] + 0.025;
        text = L"E", align = (:center, :bottom), fontsize = 23)

    text!(ax, p_star[1] + 0.025, p_star[2] + 0.025;
        text = L"S", align = (:left, :bottom), fontsize = 23)

    text!(ax, p_foot[1] + 0.018, p_foot[2] - 0.018;
        text = L"P", align = (:left, :top), fontsize = 18)

    text!(ax, O[1] - 0.018, O[2] - 0.018;
        text = L"O", align = (:right, :top), fontsize = 17)

    text!(ax, λ_mid[1], λ_mid[2] - 0.038;
        text = L"\lambda", align = (:center, :top), fontsize = 22)

    text!(ax, β_mid[1] + 0.040, β_mid[2];
        text = L"\beta", align = (:left, :center), fontsize = 22)

    text!(ax, p_zero[1] - 0.022, p_zero[2] - 0.025;
        text = L"\lambda = 0", align = (:right, :top), fontsize = 16)

    text!(ax,
        ecliptic_label_point[1] - 0.02,
        ecliptic_label_point[2] - 0.035;
        text = "Ecliptic",
        align = (:right, :top),
        fontsize = 18,
    )

    text!(ax,
        equator_label_point[1] + 0.015,
        equator_label_point[2] - 0.025;
        text = "Equator",
        align = (:left, :top),
        fontsize = 18,
    )

    return fig
end

fig = with_theme(theme_latexfonts()) do
    ecliptic_coordinate_figure()
end

pdf_path = joinpath(@__DIR__, "ecliptic_coordinate_system_main_frames.pdf")
svg_path = joinpath(@__DIR__, "ecliptic_coordinate_system_main_frames.svg")

save(pdf_path, fig)
save(svg_path, fig)

println("Saved:")
println("  ", pdf_path)
println("  ", svg_path)

fig

# Clean classical-orbital-elements illustration in CairoMakie
#
# Design goals:
#   * emphasise the physical orbit and the four angular elements
#   * keep the reference plane light
#   * omit the redundant orbital-plane guide by default
#   * place inclination locally at the ascending node
#   * keep long labels outside the central construction
#
# Install once:
#   julia -e 'using Pkg; Pkg.add(["CairoMakie", "LaTeXStrings"])'
#
# Run:
#   julia orbital_elements_makie_revised.jl
#
# Output:
#   orbital_elements_makie_revised.pdf
#   orbital_elements_makie_revised.svg

using CairoMakie
using LinearAlgebra
using LaTeXStrings

# =============================================================================
# Basic vector and camera utilities
# =============================================================================

deg(x) = π * x / 180
unit(v) = Float64.(v) / norm(v)

struct OrthoCamera
    right::Vector{Float64}
    up::Vector{Float64}
    view::Vector{Float64}  # direction from the origin towards the camera
end

function camera_frame(view; world_up = [0.0, 0.0, 1.0])
    v = unit(view)
    up0 = unit(world_up)

    if abs(dot(v, up0)) > 0.98
        up0 = [0.0, 1.0, 0.0]
    end

    right = unit(cross(up0, v))
    up = unit(cross(v, right))
    return OrthoCamera(right, up, v)
end

project_point(p, camera::OrthoCamera) =
    Point2f(dot(p, camera.right), dot(p, camera.up))

camera_depth(p, camera::OrthoCamera) = dot(p, camera.view)

function rotate_about_axis(p, axis, angle)
    a = unit(axis)
    c = cos(angle)
    s = sin(angle)
    return c .* p .+ s .* cross(a, p) .+ (1 - c) * dot(a, p) .* a
end

# =============================================================================
# Orbital geometry
# =============================================================================

"""
Return the perifocal basis vectors p̂, q̂, ŵ in the inertial frame.
The adopted sequence is R₃(Ω) R₁(i) R₃(ω).
"""
function perifocal_basis(Ω, i, ω)
    cΩ, sΩ = cos(Ω), sin(Ω)
    ci, si = cos(i), sin(i)
    cω, sω = cos(ω), sin(ω)

    Q = [
        cΩ*cω - sΩ*sω*ci   -cΩ*sω - sΩ*cω*ci    sΩ*si;
        sΩ*cω + cΩ*sω*ci   -sΩ*sω + cΩ*cω*ci   -cΩ*si;
        sω*si                cω*si                 ci
    ]

    return unit(Q[:, 1]), unit(Q[:, 2]), unit(Q[:, 3])
end

function orbit_point(true_anomaly, a, e, p̂, q̂)
    semilatus = a * (1 - e^2)
    radius = semilatus / (1 + e * cos(true_anomaly))
    return radius .* (cos(true_anomaly) .* p̂ .+ sin(true_anomaly) .* q̂)
end

function orbit_curve(a, e, p̂, q̂; samples = 1200)
    anomalies = range(0, 2π; length = samples + 1)[1:end-1]
    return [orbit_point(f, a, e, p̂, q̂) for f in anomalies]
end

function circle_in_plane(u, v, radius; samples = 720)
    θ = range(0, 2π; length = samples + 1)[1:end-1]
    return [radius .* (cos(t) .* u .+ sin(t) .* v) for t in θ]
end

function signed_angle_about_axis(a, b, axis)
    aa = unit(a)
    bb = unit(b)
    n = unit(axis)
    return atan(dot(n, cross(aa, bb)), dot(aa, bb))
end

"""Angular arc centred at `centre` and lying in the plane normal to `axis`."""
function angular_arc(
        a,
        b,
        axis,
        radius;
        centre = [0.0, 0.0, 0.0],
        samples = 160,
    )

    angle = signed_angle_about_axis(a, b, axis)
    return [
        centre .+ radius .* rotate_about_axis(unit(a), axis, t)
        for t in range(0, angle; length = samples)
    ]
end

# =============================================================================
# Occlusion and depth styling
# =============================================================================

"""
Classify an arbitrary point relative to an opaque sphere under orthographic
projection.

:occluded  sphere material lies between the camera and the point
:rear      visible, but on the far side of the scene centre
:front     visible and on the camera-facing side
"""
function visibility_state(
        p,
        camera::OrthoCamera,
        sphere_radius;
        atol = 1e-9,
    )

    x = dot(p, camera.right)
    y = dot(p, camera.up)
    z = dot(p, camera.view)
    ρ² = x^2 + y^2

    if ρ² < sphere_radius^2
        z_front = sqrt(max(0.0, sphere_radius^2 - ρ²))
        if z < z_front - atol
            return :occluded
        end
    end

    return z < 0 ? :rear : :front
end

function transition_point(a, b, state_a, camera, sphere_radius; iterations = 35)
    left = copy(a)
    right = copy(b)

    for _ in 1:iterations
        middle = 0.5 .* (left .+ right)
        if visibility_state(middle, camera, sphere_radius) == state_a
            left = middle
        else
            right = middle
        end
    end

    return 0.5 .* (left .+ right)
end

function split_curve_states(curve, camera, sphere_radius; closed = true)
    points = closed ? vcat(curve, [curve[1]]) : curve
    runs = Dict(
        :front => Vector{Vector{Point2f}}(),
        :rear => Vector{Vector{Point2f}}(),
        :occluded => Vector{Vector{Point2f}}(),
    )

    current_state = visibility_state(points[1], camera, sphere_radius)
    current_run = Point2f[project_point(points[1], camera)]

    for k in 1:(length(points) - 1)
        a = points[k]
        b = points[k + 1]
        state_a = visibility_state(a, camera, sphere_radius)
        state_b = visibility_state(b, camera, sphere_radius)

        if state_a == state_b
            push!(current_run, project_point(b, camera))
        else
            crossing = transition_point(a, b, state_a, camera, sphere_radius)
            crossing2 = project_point(crossing, camera)
            push!(current_run, crossing2)

            length(current_run) > 1 && push!(runs[current_state], current_run)

            current_state = state_b
            current_run = Point2f[crossing2, project_point(b, camera)]
        end
    end

    length(current_run) > 1 && push!(runs[current_state], current_run)
    return runs
end

function draw_runs!(ax, runs, state; linewidth, color, linestyle = :solid)
    for segment in runs[state]
        lines!(
            ax,
            segment;
            linewidth = linewidth,
            color = color,
            linestyle = linestyle,
        )
    end
end

# =============================================================================
# Arrows and labels
# =============================================================================

function add_curve_arrowhead!(
        ax,
        curve2;
        backstep = 10,
        linewidth = 1.2,
        color = :black,
        tipwidth = 10,
        tiplength = 8,
    )

    length(curve2) < 2 && return nothing
    n = length(curve2)
    p1 = curve2[max(1, n - backstep)]
    p2 = curve2[end]

    arrows2d!(
        ax,
        [p1],
        [p2];
        argmode = :endpoint,
        color = color,
        shaftwidth = linewidth,
        tipwidth = tipwidth,
        tiplength = tiplength,
    )

    return nothing
end

function draw_curved_arrow!(
        ax,
        curve3,
        camera;
        linewidth = 1.45,
        color = :black,
    )

    curve2 = [project_point(p, camera) for p in curve3]
    lines!(ax, curve2; linewidth = linewidth, color = color)
    add_curve_arrowhead!(ax, curve2; linewidth = linewidth, color = color)
    return curve2
end

function label_position_on_curve(curve2; fraction = 0.55, offset = 0.045)
    n = length(curve2)
    index = clamp(round(Int, 1 + fraction * (n - 1)), 2, n - 1)

    tangent = Vec2f(curve2[index + 1] - curve2[index - 1])
    tangent /= norm(tangent)
    normal = Vec2f(-tangent[2], tangent[1])

    # Prefer the side away from the projected origin.
    if dot(normal, Vec2f(curve2[index])) < 0
        normal = -normal
    end

    # Keep all screen-space geometry in Float32 / Point2f. Mixing Point2f
    # with Point2d creates Vector{Point2}, which recent CairoMakie cannot
    # convert into a concrete line buffer.
    return Point2f(curve2[index] + Float32(offset) .* normal)
end

function outward_label_position(
        anchor::Point2f;
        radial = 0.14,
        tangential = 0.04,
    )

    r = Vec2f(anchor)
    if norm(r) < 1e-9
        return Point2f(anchor + Vec2f(Float32(radial), Float32(tangential)))
    end

    r /= norm(r)
    t = Vec2f(-r[2], r[1])
    return Point2f(
        anchor + Float32(radial) .* r + Float32(tangential) .* t
    )
end

function leader_label!(
        ax,
        anchor::Point2f,
        label_position::Point2f,
        text_value;
        align = (:left, :center),
        fontsize = 15,
        elbow_fraction = 0.55,
    )

    delta = Vec2f(label_position - anchor)
    elbow = Point2f(anchor + Float32(elbow_fraction) .* delta)

    # A short orthogonal bend usually avoids the appearance of one long ray.
    tangent = norm(delta) > 1f-9 ? delta / norm(delta) : Vec2f(1f0, 0f0)
    normal = Vec2f(-tangent[2], tangent[1])
    elbow = Point2f(elbow + 0.018f0 .* normal)

    # Explicit concrete element type avoids a heterogeneous Vector{Point2}.
    leader_path = Point2f[anchor, elbow, label_position]
    lines!(
        ax,
        leader_path;
        color = (:black, 0.48),
        linewidth = 0.65,
    )

    text!(
        ax,
        label_position[1],
        label_position[2];
        text = text_value,
        align = align,
        fontsize = fontsize,
    )
end

# =============================================================================
# Main figure
# =============================================================================

function orbital_elements_figure(;
        a = 1.0,
        e = 0.32,
        Ω = deg(38),
        i = deg(34),
        ω = deg(52),
        ν = deg(78),
        earth_radius = 0.16,
        camera_view = [-3.4, -4.5, 2.9],
        show_equatorial_guide = true,
        show_orbital_plane_guide = false,
        show_line_of_nodes_label = false,
        show_orbit_label = true,
    )

    camera = camera_frame(camera_view)

    # Inertial reference frame.
    x̂ = [1.0, 0.0, 0.0]
    ŷ = [0.0, 1.0, 0.0]
    ẑ = [0.0, 0.0, 1.0]

    # Orbital/perifocal frame.
    p̂, q̂, ŵ = perifocal_basis(Ω, i, ω)
    n̂ = unit([cos(Ω), sin(Ω), 0.0])
    t̂eq = unit(cross(ẑ, n̂))
    t̂orb = unit(cross(ŵ, n̂))

    orbit = orbit_curve(a, e, p̂, q̂)
    satellite = orbit_point(ν, a, e, p̂, q̂)
    perigee = orbit_point(0.0, a, e, p̂, q̂)
    ascending_node = orbit_point(-ω, a, e, p̂, q̂)

    guide_radius = 1.30 * a
    equatorial_guide = circle_in_plane(x̂, ŷ, guide_radius)
    orbital_guide = circle_in_plane(n̂, t̂orb, guide_radius)

    orbit_runs = split_curve_states(orbit, camera, earth_radius)
    eq_runs = split_curve_states(equatorial_guide, camera, earth_radius)
    orbplane_runs = split_curve_states(orbital_guide, camera, earth_radius)

    # Distinct radii keep the angular annotations visually separated.
    Ω_arc = angular_arc(x̂, n̂, ẑ, 0.46)
    ω_arc = angular_arc(n̂, p̂, ŵ, 0.63)
    ν_arc = angular_arc(p̂, unit(satellite), ŵ, 0.82)

    # Inclination is shown locally at the ascending node instead of at Earth.
    i_centre = 0.74 .* ascending_node
    i_arc = angular_arc(
        t̂eq,
        t̂orb,
        n̂,
        0.19;
        centre = i_centre,
    )

    fig = Figure(size = (980, 790), figure_padding = 12)
    ax = Axis(fig[1, 1]; aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    xlims!(ax, -1.58, 1.58)
    ylims!(ax, -1.27, 1.27)

    # -------------------------------------------------------------------------
    # Rear geometry
    # -------------------------------------------------------------------------

    if show_equatorial_guide
        draw_runs!(
            ax, eq_runs, :rear;
            linewidth = 0.72,
            color = (:black, 0.17),
            linestyle = (:dash, :loose),
        )
    end

    if show_orbital_plane_guide
        draw_runs!(
            ax, orbplane_runs, :rear;
            linewidth = 0.62,
            color = (:black, 0.13),
            linestyle = (:dash, :loose),
        )
    end

    draw_runs!(
        ax, orbit_runs, :rear;
        linewidth = 1.55,
        color = (:black, 0.42),
        linestyle = (:dash, :loose),
    )

    # Line of nodes: present, but not visually dominant.
    node_start = project_point(-1.38 * a .* n̂, camera)
    node_end = project_point(1.38 * a .* n̂, camera)
    lines!(ax, Point2f[node_start, node_end]; color = (:black, 0.34), linewidth = 0.85)

    # -------------------------------------------------------------------------
    # Opaque Earth
    # -------------------------------------------------------------------------

    earth_outline = [
        Point2f(earth_radius * cos(t), earth_radius * sin(t))
        for t in range(0, 2π; length = 360)
    ]

    poly!(
        ax,
        earth_outline;
        color = :white,
        strokecolor = :black,
        strokewidth = 1.45,
    )

    earth_equator = [
        Point2f(earth_radius * cos(t), 0.28 * earth_radius * sin(t))
        for t in range(0, 2π; length = 240)
    ]
    lines!(ax, earth_equator; color = (:black, 0.30), linewidth = 0.65)

    # -------------------------------------------------------------------------
    # Front geometry
    # -------------------------------------------------------------------------

    if show_equatorial_guide
        draw_runs!(
            ax, eq_runs, :front;
            linewidth = 0.82,
            color = (:black, 0.32),
        )
    end

    if show_orbital_plane_guide
        draw_runs!(
            ax, orbplane_runs, :front;
            linewidth = 0.72,
            color = (:black, 0.22),
        )
    end

    draw_runs!(
        ax, orbit_runs, :front;
        linewidth = 2.05,
        color = :black,
    )

    O = Point2f(0, 0)

    # Reference direction.
    ref_end = project_point(1.42 * a .* x̂, camera)
    arrows2d!(
        ax, [O], [ref_end];
        argmode = :endpoint,
        color = (:black, 0.58),
        shaftwidth = 0.9,
        tipwidth = 9,
        tiplength = 7,
    )

    # Satellite radius vector.
    p_sat = project_point(satellite, camera)
    arrows2d!(
        ax, [O], [p_sat];
        argmode = :endpoint,
        color = :black,
        shaftwidth = 1.15,
        tipwidth = 10,
        tiplength = 8,
    )

    # Angular arrows are explanatory overlays and intentionally remain visible.
    Ω2 = draw_curved_arrow!(ax, Ω_arc, camera; linewidth = 1.35)
    i2 = draw_curved_arrow!(ax, i_arc, camera; linewidth = 1.30)
    ω2 = draw_curved_arrow!(ax, ω_arc, camera; linewidth = 1.35)
    ν2 = draw_curved_arrow!(ax, ν_arc, camera; linewidth = 1.35)

    p_perigee = project_point(perigee, camera)
    p_node = project_point(ascending_node, camera)

    scatter!(ax, [O]; color = :black, markersize = 5)
    scatter!(ax, [p_node]; color = :black, markersize = 7)
    scatter!(ax, [p_perigee]; color = :black, markersize = 7)
    scatter!(ax, [p_sat]; color = :black, markersize = 9)

    # -------------------------------------------------------------------------
    # Mathematical labels
    # -------------------------------------------------------------------------

    Ω_label = label_position_on_curve(Ω2; fraction = 0.52, offset = 0.050)
    i_label = label_position_on_curve(i2; fraction = 0.52, offset = 0.040)
    ω_label = label_position_on_curve(ω2; fraction = 0.54, offset = 0.055)
    ν_label = label_position_on_curve(ν2; fraction = 0.58, offset = 0.060)

    # Named screen-space corrections are easy to tune after changing the camera.
    Ω_label += Vec2f(0.000, 0.012)
    i_label += Vec2f(-0.010, 0.000)
    ω_label += Vec2f(0.015, 0.012)
    ν_label += Vec2f(-0.012, 0.010)

    text!(ax, Ω_label[1], Ω_label[2]; text = L"\Omega", fontsize = 22,
        align = (:center, :center))
    text!(ax, i_label[1], i_label[2]; text = L"i", fontsize = 22,
        align = (:center, :center))
    text!(ax, ω_label[1], ω_label[2]; text = L"\omega", fontsize = 22,
        align = (:center, :center))
    text!(ax, ν_label[1], ν_label[2]; text = L"\nu", fontsize = 22,
        align = (:center, :center))

    text!(ax, 0, 0; text = "Earth", fontsize = 12,
        align = (:center, :center))

    # -------------------------------------------------------------------------
    # Restrained descriptive labels
    # -------------------------------------------------------------------------

    satellite_label = outward_label_position(p_sat; radial = 0.13, tangential = 0.05)
    perigee_label = outward_label_position(p_perigee; radial = 0.14, tangential = -0.04)
    node_label = outward_label_position(p_node; radial = 0.13, tangential = -0.05)

    leader_label!(
        ax,
        p_sat,
        satellite_label,
        "Satellite";
        align = (:left, :bottom),
        fontsize = 15,
    )

    leader_label!(
        ax,
        p_perigee,
        perigee_label,
        "Perigee";
        align = (:left, :bottom),
        fontsize = 15,
    )

    leader_label!(
        ax,
        p_node,
        node_label,
        "Ascending node";
        align = (:left, :top),
        fontsize = 15,
    )

    if show_line_of_nodes_label
        line_nodes_anchor = project_point(1.08 * a .* n̂, camera)
        line_nodes_label = outward_label_position(
            line_nodes_anchor;
            radial = 0.11,
            tangential = 0.05,
        )
        leader_label!(
            ax,
            line_nodes_anchor,
            line_nodes_label,
            "Line of nodes";
            align = (:left, :bottom),
            fontsize = 14,
        )
    end

    text!(
        ax,
        ref_end[1] + 0.025,
        ref_end[2] - 0.025;
        text = "Reference direction",
        fontsize = 14,
        align = (:left, :top),
    )

    if show_equatorial_guide
        eq_label_anchor = project_point(
            guide_radius .* (cos(deg(155)) .* x̂ .+ sin(deg(155)) .* ŷ),
            camera,
        )
        text!(
            ax,
            eq_label_anchor[1],
            eq_label_anchor[2] - 0.030;
            text = "Equatorial plane",
            fontsize = 14,
            color = (:black, 0.68),
            align = (:center, :top),
        )
    end

    if show_orbital_plane_guide
        orbplane_label_anchor = project_point(
            guide_radius .* (cos(deg(232)) .* n̂ .+ sin(deg(232)) .* t̂orb),
            camera,
        )
        text!(
            ax,
            orbplane_label_anchor[1],
            orbplane_label_anchor[2] - 0.030;
            text = "Orbital plane",
            fontsize = 14,
            color = (:black, 0.55),
            align = (:center, :top),
        )
    end

    if show_orbit_label
        orbit_label_anchor = project_point(
            orbit_point(deg(220), a, e, p̂, q̂),
            camera,
        )
        text!(
            ax,
            orbit_label_anchor[1],
            orbit_label_anchor[2] - 0.045;
            text = "Satellite orbit",
            fontsize = 14,
            align = (:center, :top),
        )
    end

    return fig
end

fig = with_theme(theme_latexfonts()) do
    orbital_elements_figure()
end

pdf_path = joinpath(@__DIR__, "orbital_elements_makie_revised.pdf")
svg_path = joinpath(@__DIR__, "orbital_elements_makie_revised.svg")

save(pdf_path, fig)
save(svg_path, fig)

println("Saved:")
println("  ", pdf_path)
println("  ", svg_path)

fig

# Composed classical-orbital-elements illustration in CairoMakie
#
# This version revises the scene composition:
#   * the reference direction is brought to the visual foreground
#   * perigee is placed on the front/upper side of the orbit and labelled outside
#   * the inclination arc is shown at the front node for clarity
#   * the satellite is a larger square marker
#   * Earth and arrows are drawn more cleanly
#
# Install once:
#   julia -e 'using Pkg; Pkg.add(["CairoMakie", "LaTeXStrings"])'
#
# Run:
#   julia orbital_elements_makie_composed.jl
#
# Output:
#   orbital_elements_makie_composed.pdf
#   orbital_elements_makie_composed.svg

using CairoMakie
using LinearAlgebra
using LaTeXStrings

# =============================================================================
# Basic vector and camera utilities
# =============================================================================

deg(x) = π * x / 180
unit(v) = Float64.(v) / norm(v)

struct OrthoCamera
    right::Vector{Float64}
    up::Vector{Float64}
    view::Vector{Float64}  # direction from the origin towards the camera
end

function camera_frame(view; world_up = [0.0, 0.0, 1.0])
    v = unit(view)
    up0 = unit(world_up)

    if abs(dot(v, up0)) > 0.98
        up0 = [0.0, 1.0, 0.0]
    end

    right = unit(cross(up0, v))
    up = unit(cross(v, right))
    return OrthoCamera(right, up, v)
end

project_point(p, camera::OrthoCamera) =
    Point2f(dot(p, camera.right), dot(p, camera.up))

camera_depth(p, camera::OrthoCamera) = dot(p, camera.view)

function rotate_about_axis(p, axis, angle)
    a = unit(axis)
    c = cos(angle)
    s = sin(angle)
    return c .* p .+ s .* cross(a, p) .+ (1 - c) * dot(a, p) .* a
end

# =============================================================================
# Orbital geometry
# =============================================================================

function perifocal_basis(Ω, i, ω)
    cΩ, sΩ = cos(Ω), sin(Ω)
    ci, si = cos(i), sin(i)
    cω, sω = cos(ω), sin(ω)

    Q = [
        cΩ*cω - sΩ*sω*ci   -cΩ*sω - sΩ*cω*ci    sΩ*si;
        sΩ*cω + cΩ*sω*ci   -sΩ*sω + cΩ*cω*ci   -cΩ*si;
        sω*si                cω*si                 ci
    ]

    return unit(Q[:, 1]), unit(Q[:, 2]), unit(Q[:, 3])
end

function orbit_point(true_anomaly, a, e, p̂, q̂)
    semilatus = a * (1 - e^2)
    radius = semilatus / (1 + e * cos(true_anomaly))
    return radius .* (cos(true_anomaly) .* p̂ .+ sin(true_anomaly) .* q̂)
end

function orbit_curve(a, e, p̂, q̂; samples = 1200)
    anomalies = range(0, 2π; length = samples + 1)[1:end-1]
    return [orbit_point(f, a, e, p̂, q̂) for f in anomalies]
end

function circle_in_plane(u, v, radius; samples = 720)
    θ = range(0, 2π; length = samples + 1)[1:end-1]
    return [radius .* (cos(t) .* u .+ sin(t) .* v) for t in θ]
end

function signed_angle_about_axis(a, b, axis)
    aa = unit(a)
    bb = unit(b)
    n = unit(axis)
    return atan(dot(n, cross(aa, bb)), dot(aa, bb))
end

function angular_arc(
        a,
        b,
        axis,
        radius;
        centre = [0.0, 0.0, 0.0],
        samples = 160,
    )

    angle = signed_angle_about_axis(a, b, axis)
    return [
        centre .+ radius .* rotate_about_axis(unit(a), axis, t)
        for t in range(0, angle; length = samples)
    ]
end

# =============================================================================
# Occlusion and depth styling
# =============================================================================

function visibility_state(
        p,
        camera::OrthoCamera,
        sphere_radius;
        atol = 1e-9,
    )

    x = dot(p, camera.right)
    y = dot(p, camera.up)
    z = dot(p, camera.view)
    ρ² = x^2 + y^2

    if ρ² < sphere_radius^2
        z_front = sqrt(max(0.0, sphere_radius^2 - ρ²))
        if z < z_front - atol
            return :occluded
        end
    end

    return z < 0 ? :rear : :front
end

function transition_point(a, b, state_a, camera, sphere_radius; iterations = 35)
    left = copy(a)
    right = copy(b)

    for _ in 1:iterations
        middle = 0.5 .* (left .+ right)
        if visibility_state(middle, camera, sphere_radius) == state_a
            left = middle
        else
            right = middle
        end
    end

    return 0.5 .* (left .+ right)
end

function split_curve_states(curve, camera, sphere_radius; closed = true)
    points = closed ? vcat(curve, [curve[1]]) : curve
    runs = Dict(
        :front => Vector{Vector{Point2f}}(),
        :rear => Vector{Vector{Point2f}}(),
        :occluded => Vector{Vector{Point2f}}(),
    )

    current_state = visibility_state(points[1], camera, sphere_radius)
    current_run = Point2f[project_point(points[1], camera)]

    for k in 1:(length(points) - 1)
        a = points[k]
        b = points[k + 1]
        state_a = visibility_state(a, camera, sphere_radius)
        state_b = visibility_state(b, camera, sphere_radius)

        if state_a == state_b
            push!(current_run, project_point(b, camera))
        else
            crossing = transition_point(a, b, state_a, camera, sphere_radius)
            crossing2 = project_point(crossing, camera)
            push!(current_run, crossing2)
            length(current_run) > 1 && push!(runs[current_state], current_run)
            current_state = state_b
            current_run = Point2f[crossing2, project_point(b, camera)]
        end
    end

    length(current_run) > 1 && push!(runs[current_state], current_run)
    return runs
end

function draw_runs!(ax, runs, state; linewidth, color, linestyle = :solid)
    for segment in runs[state]
        lines!(ax, segment; linewidth = linewidth, color = color, linestyle = linestyle)
    end
end

# =============================================================================
# Arrows and labels
# =============================================================================

function add_curve_arrowhead!(
        ax,
        curve2;
        backstep = 10,
        linewidth = 1.2,
        color = :black,
        tipwidth = 10,
        tiplength = 8,
    )

    length(curve2) < 2 && return nothing
    n = length(curve2)
    p1 = curve2[max(1, n - backstep)]
    p2 = curve2[end]

    arrows2d!(
        ax,
        [p1],
        [p2];
        argmode = :endpoint,
        color = color,
        shaftwidth = linewidth,
        tipwidth = tipwidth,
        tiplength = tiplength,
    )
    return nothing
end

function draw_curved_arrow!(ax, curve3, camera; linewidth = 1.45, color = :black)
    curve2 = [project_point(p, camera) for p in curve3]
    lines!(ax, curve2; linewidth = linewidth, color = color)
    add_curve_arrowhead!(ax, curve2; linewidth = linewidth, color = color)
    return curve2
end

function label_position_on_curve(curve2; fraction = 0.55, offset = 0.045)
    n = length(curve2)
    index = clamp(round(Int, 1 + fraction * (n - 1)), 2, n - 1)
    tangent = Vec2f(curve2[index + 1] - curve2[index - 1])
    tangent /= norm(tangent)
    normal = Vec2f(-tangent[2], tangent[1])
    if dot(normal, Vec2f(curve2[index])) < 0
        normal = -normal
    end
    return Point2f(curve2[index] + Float32(offset) .* normal)
end

function outward_label_position(anchor::Point2f; radial = 0.14, tangential = 0.04)
    r = Vec2f(anchor)
    if norm(r) < 1e-9
        return Point2f(anchor + Vec2f(Float32(radial), Float32(tangential)))
    end
    r /= norm(r)
    t = Vec2f(-r[2], r[1])
    return Point2f(anchor + Float32(radial) .* r + Float32(tangential) .* t)
end

function leader_label!(
        ax,
        anchor::Point2f,
        label_position::Point2f,
        text_value;
        align = (:left, :center),
        fontsize = 15,
        elbow_fraction = 0.55,
        color = (:black, 0.50),
        linewidth = 0.65,
    )

    delta = Vec2f(label_position - anchor)
    elbow = Point2f(anchor + Float32(elbow_fraction) .* delta)
    tangent = norm(delta) > 1f-9 ? delta / norm(delta) : Vec2f(1f0, 0f0)
    normal = Vec2f(-tangent[2], tangent[1])
    elbow = Point2f(elbow + 0.016f0 .* normal)

    leader_path = Point2f[anchor, elbow, label_position]
    lines!(ax, leader_path; color = color, linewidth = linewidth)
    text!(ax, label_position[1], label_position[2];
        text = text_value, align = align, fontsize = fontsize)
end

# =============================================================================
# Main figure
# =============================================================================

function orbital_elements_figure(;
        a = 1.0,
        e = 0.38,
        Ω = deg(34),
        i = deg(42),
        ω = deg(110),
        ν = deg(52),
        earth_radius = 0.17,
        camera_view = [3.8, -2.6, 2.6],
        show_equatorial_guide = true,
        show_orbital_plane_guide = false,
        show_line_of_nodes_label = false,
        show_orbit_label = true,
    )

    camera = camera_frame(camera_view)

    x̂ = [1.0, 0.0, 0.0]
    ŷ = [0.0, 1.0, 0.0]
    ẑ = [0.0, 0.0, 1.0]

    p̂, q̂, ŵ = perifocal_basis(Ω, i, ω)
    n̂ = unit([cos(Ω), sin(Ω), 0.0])
    t̂eq = unit(cross(ẑ, n̂))
    t̂orb = unit(cross(ŵ, n̂))

    orbit = orbit_curve(a, e, p̂, q̂)
    satellite = orbit_point(ν, a, e, p̂, q̂)
    perigee = orbit_point(0.0, a, e, p̂, q̂)
    ascending_node = orbit_point(-ω, a, e, p̂, q̂)

    # Use the node that is visually closer to the viewer for the inclination arc.
    descending_node = -ascending_node
    i_node = camera_depth(ascending_node, camera) >= camera_depth(descending_node, camera) ? ascending_node : descending_node

    guide_radius = 1.30 * a
    equatorial_guide = circle_in_plane(x̂, ŷ, guide_radius)
    orbital_guide = circle_in_plane(n̂, t̂orb, guide_radius)

    orbit_runs = split_curve_states(orbit, camera, earth_radius)
    eq_runs = split_curve_states(equatorial_guide, camera, earth_radius)
    orbplane_runs = split_curve_states(orbital_guide, camera, earth_radius)

    # Distinct radii keep the angular annotations visually separated.
    Ω_arc = angular_arc(x̂, n̂, ẑ, 0.50)
    ω_arc = angular_arc(n̂, p̂, ŵ, 0.66)
    ν_arc = angular_arc(p̂, unit(satellite), ŵ, 0.88)

    # Inclination shown on the visually front node.
    i_centre = 0.82 .* i_node
    i_arc = angular_arc(t̂eq, t̂orb, n̂, 0.22; centre = i_centre)

    fig = Figure(size = (1020, 800), figure_padding = 12)
    ax = Axis(fig[1, 1]; aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    xlims!(ax, -1.62, 1.62)
    ylims!(ax, -1.20, 1.28)

    # -------------------------------------------------------------------------
    # Rear geometry
    # -------------------------------------------------------------------------
    if show_equatorial_guide
        draw_runs!(ax, eq_runs, :rear;
            linewidth = 0.78,
            color = (:black, 0.16),
            linestyle = (:dash, :loose),
        )
    end

    if show_orbital_plane_guide
        draw_runs!(ax, orbplane_runs, :rear;
            linewidth = 0.62,
            color = (:black, 0.12),
            linestyle = (:dash, :loose),
        )
    end

    draw_runs!(ax, orbit_runs, :rear;
        linewidth = 1.55,
        color = (:black, 0.42),
        linestyle = (:dash, :loose),
    )

    node_start = project_point(-1.35 * a .* n̂, camera)
    node_end = project_point(1.35 * a .* n̂, camera)
    lines!(ax, Point2f[node_start, node_end]; color = (:black, 0.28), linewidth = 0.82)

    # -------------------------------------------------------------------------
    # Opaque Earth
    # -------------------------------------------------------------------------
    earth_outline = [Point2f(earth_radius * cos(t), earth_radius * sin(t)) for t in range(0, 2π; length = 360)]
    poly!(ax, earth_outline; color = RGBf(0.98, 0.98, 0.98), strokecolor = :black, strokewidth = 1.45)

    earth_equator = [Point2f(earth_radius * cos(t), 0.30 * earth_radius * sin(t)) for t in range(0, 2π; length = 240)]
    lines!(ax, earth_equator; color = (:black, 0.28), linewidth = 0.68)

    # -------------------------------------------------------------------------
    # Front geometry
    # -------------------------------------------------------------------------
    if show_equatorial_guide
        draw_runs!(ax, eq_runs, :front;
            linewidth = 0.92,
            color = (:black, 0.34),
        )
    end

    if show_orbital_plane_guide
        draw_runs!(ax, orbplane_runs, :front;
            linewidth = 0.72,
            color = (:black, 0.22),
        )
    end

    draw_runs!(ax, orbit_runs, :front; linewidth = 2.10, color = :black)

    O = Point2f(0, 0)

    # Reference direction intentionally in the visual foreground.
    ref_end = project_point(1.55 * a .* x̂, camera)
    arrows2d!(ax, [O], [ref_end];
        argmode = :endpoint,
        color = (:black, 0.72),
        shaftwidth = 1.0,
        tipwidth = 10,
        tiplength = 8,
    )

    # Satellite radius vector.
    p_sat = project_point(satellite, camera)
    arrows2d!(ax, [O], [p_sat];
        argmode = :endpoint,
        color = :black,
        shaftwidth = 1.20,
        tipwidth = 10,
        tiplength = 8,
    )

    # Angular arrows.
    Ω2 = draw_curved_arrow!(ax, Ω_arc, camera; linewidth = 1.38)
    i2 = draw_curved_arrow!(ax, i_arc, camera; linewidth = 1.34)
    ω2 = draw_curved_arrow!(ax, ω_arc, camera; linewidth = 1.38)
    ν2 = draw_curved_arrow!(ax, ν_arc, camera; linewidth = 1.38)

    p_perigee = project_point(perigee, camera)
    p_node = project_point(ascending_node, camera)

    scatter!(ax, [O]; color = :black, markersize = 5)
    scatter!(ax, [p_node]; color = :black, markersize = 7)
    scatter!(ax, [p_perigee]; color = :black, markersize = 7)
    scatter!(ax, [p_sat]; color = :black, markersize = 16, marker = :rect)

    # -------------------------------------------------------------------------
    # Mathematical labels
    # -------------------------------------------------------------------------
    Ω_label = label_position_on_curve(Ω2; fraction = 0.54, offset = 0.050)
    i_label = label_position_on_curve(i2; fraction = 0.50, offset = 0.038)
    ω_label = label_position_on_curve(ω2; fraction = 0.56, offset = 0.055)
    ν_label = label_position_on_curve(ν2; fraction = 0.40, offset = 0.060)

    Ω_label += Vec2f(0.000, 0.012)
    i_label += Vec2f(-0.006, 0.000)
    ω_label += Vec2f(0.016, 0.012)
    ν_label += Vec2f(-0.020, 0.006)

    text!(ax, Ω_label[1], Ω_label[2]; text = L"\Omega", fontsize = 22, align = (:center, :center))
    text!(ax, i_label[1], i_label[2]; text = L"i", fontsize = 22, align = (:center, :center))
    text!(ax, ω_label[1], ω_label[2]; text = L"\omega", fontsize = 22, align = (:center, :center))
    text!(ax, ν_label[1], ν_label[2]; text = L"\nu", fontsize = 22, align = (:center, :center))

    text!(ax, 0, 0; text = "Earth", fontsize = 12, align = (:center, :center))

    # -------------------------------------------------------------------------
    # Restrained descriptive labels
    # -------------------------------------------------------------------------
    satellite_label = outward_label_position(p_sat; radial = 0.14, tangential = 0.05)
    # Perigee label intentionally moved outside / above the guide plane region.
    perigee_label = Point2f(p_perigee + Vec2f(0.02f0, 0.12f0))
    node_label = outward_label_position(p_node; radial = 0.12, tangential = -0.06)

    leader_label!(ax, p_sat, satellite_label, "Satellite";
        align = (:right, :center), fontsize = 15, elbow_fraction = 0.52)

    leader_label!(ax, p_perigee, perigee_label, "Perigee";
        align = (:left, :bottom), fontsize = 15, elbow_fraction = 0.42)

    leader_label!(ax, p_node, node_label, "Ascending node";
        align = (:left, :center), fontsize = 15, elbow_fraction = 0.54)

    if show_line_of_nodes_label
        line_nodes_anchor = project_point(1.02 * a .* n̂, camera)
        line_nodes_label = outward_label_position(line_nodes_anchor; radial = 0.10, tangential = 0.04)
        leader_label!(ax, line_nodes_anchor, line_nodes_label, "Line of nodes";
            align = (:left, :bottom), fontsize = 14)
    end

    text!(ax, ref_end[1] + 0.025, ref_end[2] - 0.020;
        text = "Reference direction", fontsize = 14, align = (:left, :top))

    if show_equatorial_guide
        eq_label_anchor = project_point(guide_radius .* (cos(deg(170)) .* x̂ .+ sin(deg(170)) .* ŷ), camera)
        text!(ax, eq_label_anchor[1], eq_label_anchor[2] - 0.030;
            text = "Equatorial plane", fontsize = 14, color = (:black, 0.66), align = (:center, :top))
    end

    if show_orbital_plane_guide
        orbplane_label_anchor = project_point(guide_radius .* (cos(deg(232)) .* n̂ .+ sin(deg(232)) .* t̂orb), camera)
        text!(ax, orbplane_label_anchor[1], orbplane_label_anchor[2] - 0.030;
            text = "Orbital plane", fontsize = 14, color = (:black, 0.55), align = (:center, :top))
    end

    if show_orbit_label
        orbit_label_anchor = project_point(orbit_point(deg(230), a, e, p̂, q̂), camera)
        text!(ax, orbit_label_anchor[1], orbit_label_anchor[2] - 0.050;
            text = "Satellite orbit", fontsize = 14, align = (:center, :top))
    end

    return fig
end

fig = with_theme(theme_latexfonts()) do
    orbital_elements_figure()
end

pdf_path = joinpath(@__DIR__, "orbital_elements_makie_composed.pdf")
svg_path = joinpath(@__DIR__, "orbital_elements_makie_composed.svg")

save(pdf_path, fig)
save(svg_path, fig)

println("Saved:")
println("  ", pdf_path)
println("  ", svg_path)

fig