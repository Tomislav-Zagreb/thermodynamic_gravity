import numpy as np
import matplotlib.pyplot as plt

########################################################
########################################################

def velocity_profile(radiusList, totalMass, systemRadius, alpha=0.2, n=0.01):
    """
    Calculate the velocity profile of a galaxy based on given parametres.

    Parametres:
    radiusList (array): Radial distances (in parsecs) at which the velocity is calculated.
    totalMass (float): Total mass of the galaxy (in solar masses).
    systemRadius (float): Total radius of the system like a galaxy (in parsecs).
    alpha (float): Position of the peak velocity (fraction of d, typically between 0 and 1).
    n (float): Decay parametre for the velocity drop after the peak.

    Returns:
    velocities (array): Corresponding velocities at distances r.
    """
    # Gravitational constant (in units compatible with solar mass and parsecs)
    # [pc * M_sun^-1 * (km/s)^2]
    G = 4.30091e-3 

    # Maximum velocity at r = alpha * d
    vmax = np.sqrt(G * M / (alpha * d))

    # Calculate velocity profile
    velocities = np.zeros_like(r)

    for i, ri in enumerate(r):
        if ri <= alpha * d:
            # Inside the galaxy, increasing velocity
            velocities[i] = vmax * (ri / (alpha * d)) * np.exp(-ri / (alpha * d))
        else:
            # Outside the galaxy, decreasing velocity
            v_at_alpha_d = vmax * np.exp(-1)  # Value of velocity at r = alpha * d
            velocities[i] = v_at_alpha_d * np.exp(-n * (ri / d - alpha))

    return velocities

########################################################
########################################################

M = 1e11    # Total mass of the galaxy in solar masses 
d = 30000   # Total radius of the galaxy in parsecs  
n = 0.1     # Decay parametre for the flat decrease
# alpha = 0.2  # Location of the peak velocity (20 % of d)

# Generate radial distances for the plot (from 0 to factor * d)
r = np.linspace(1, 2.00 * d, 1000)

########################################################
########################################################

def calculate_alpha(M, d, K_alpha=0.03, massScalingFactor=1e6, distanceScalingFactor=1e3): 
    """
    Calculate alpha (peak velocity radius as a fraction of d).

    Parametres:
    M : float : Mass of the galaxy (solar masses)
    d : float : Diameter of the galaxy (parsecs)
    K_alpha : float : Scaling constant

    Returns:
    float : Alpha value clamped between 0 and 0.3
    """
    # Compute alpha using logarithmic scaling
    alpha = K_alpha * np.log10((M / massScalingFactor) * (d / distanceScalingFactor))

    # Clamp alpha between 0 and 0.3
    return max(0.0001, min(0.3, alpha))

alphaNew = calculate_alpha(M, d)

########################################################
########################################################

def calculate_n(M, d, massScalingFactor=1e6, distanceScalingFactor=1e3, s=7.75):
    """
    Calculate n (flat decrease decay parametre) based on galaxy mass and diametre.

    Parametres:
    M : float : Mass of the galaxy (solar masses)
    d : float : Diameter of the galaxy (parsecs)
    M0 : float : Reference mass (e.g., 1e6 solar masses)
    d0 : float : Reference diameter (e.g., 1e3 parsecs)
    s : float : Scaling factor for the rate of change

    Returns:
    float : n value clamped between 0 and 1
    """
    # Calculate n
    n = 1 - (np.log10(M / massScalingFactor) + np.log10(d / distanceScalingFactor)) / s
    # Clamp n between 0 and 1
    return max(0.01, min(1, n))

nNew = calculate_n(M, d)

########################################################
########################################################
# Compute the velocity profile

cosmicScalingCoefficientForFinalVelocities = 2.5
v = velocity_profile(r, M, d, alphaNew, nNew) * cosmicScalingCoefficientForFinalVelocities


########################################################
########################################################
# Plot the velocity profile

plt.figure(figsize=(8, 6))
plt.plot(r, v, label=f'Alpha = {alphaNew}', color='blue')
plt.axvline(x=alphaNew*d, color='red', linestyle='--', label=f'Peak at {alphaNew*d:.0f} pc, n is {nNew}')
plt.xlabel('Radius (pc)')
plt.ylabel('Velocity (km/s)')
plt.title('Galaxy Velocity Profile')
plt.legend()
plt.grid(True)
plt.show()

########################################################
########################################################
