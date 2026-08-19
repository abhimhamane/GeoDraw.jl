# keplerian elements

circle_point(longitude; radius = 1.0) =
    radius .* [cos(longitude), sin(longitude), 0.0]

# ECI frame
x_eci = unit(circle_point(120.0))
z_eci = [0.0, 0.0, 1.0]
y_eci = cross(z_eci, x_eci)

