function main_0D(du, u, p, t)

    b = u[1]
    w = u[2]

    du[1] = p[1] * b * (1-b) * ((1 + p[2]*b)^2) * w - b
    du[2] = p[3] - p[1]*w - p[4]*((1 + p[2]*b)^2) * b * w

end
