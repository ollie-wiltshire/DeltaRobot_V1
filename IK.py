import math


def delta_ik(
    x: float,
    y: float,
    z: float,
    r_b: float,
    r_p: float,
    L1: float,
    L2: float,
    phi_deg: tuple[float, float, float] = (0.0, 120.0, 240.0),
    branch: int = 1,
) -> tuple[float, float, float]:
    """
    Compute joint angles (theta_1, theta_2, theta_3) in radians for a 3-DOF Delta robot
    given target platform position (x, y, z) with +z pointing downward.
    """
   
    if branch not in (1, -1):
        raise ValueError("branch must be +1 or -1")


    d = r_b - r_p
    thetas: list[float] = []

    for phi_d in phi_deg:
        phi = math.radians(phi_d)
        c_phi = math.cos(phi)
        s_phi = math.sin(phi)

        # Rotate global (x, y, z) into arm i's local (u, v, w) frame
        u_i = x * c_phi + y * s_phi
        v_i = -x * s_phi + y * c_phi

        du_i = u_i - d

        # Coefficients for A_i * cos(theta_i) + B_i * sin(theta_i) = C_i
        A_i = -2.0 * L1 * du_i
        B_i = 2.0 * L1 * z
        C_i = L2**2 - (du_i**2 + v_i**2 + z**2 + L1**2)

        discriminant = A_i**2 + B_i**2 - C_i**2
        if discriminant < 0.0:
            raise ValueError(
                f"Target ({x}, {y}, {z}) is unreachable for arm at phi={phi_d} deg "
                f"(discriminant={discriminant:.4f})"
            )

        sqrt_disc = math.sqrt(discriminant)
        # t = tan(theta_i / 2) = (B_i +/- sqrt(A_i^2 + B_i^2 - C_i^2)) / (C_i + A_i)
        # theta_i = 2 * atan2(numerator, denominator)
        numerator = B_i + branch * sqrt_disc
        denominator = C_i + A_i
        theta_i = 2.0 * math.atan2(numerator, denominator)

        thetas.append(theta_i)

    return (thetas[0], thetas[1], thetas[2])


if __name__ == "__main__":
    r_b = 150.0   # Base radial offset (mm)
    r_p = 50.0    # Platform radial offset (mm)
    L1 = 120.0    # Upper arm length (mm)
    L2 = 250.0    # Forearm / parallelogram length (mm)

    t1, t2, t3 = delta_ik(x=0, y=0, z=220.0, r_b=r_b, r_p=r_p, L1=L1, L2=L2)
    print("Radians:", t1, t2, t3)
    print("Degrees:", math.degrees(t1), math.degrees(t2), math.degrees(t3))