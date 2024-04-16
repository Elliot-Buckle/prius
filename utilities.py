import numpy as np
import constants

def performance(contour_mach, area_exit, velocity_exit, pressure_exit, pressure_ambient, mass_flowrate):
    """
    Calculates performance parameters along the contour
    """

    contour_temperature = temp_from_mach(contour_mach)
    contour_pressure = pressure_from_mach(contour_mach)
    contour_density = density_from_mach(
        contour_pressure, contour_temperature
    )
    contour_mdot = mass_flux(contour_mach)
    contour_velocity = velocity_from_mach(
        temperature=contour_temperature, mach=contour_mach
    )
    specific_impulse = (
        velocity_exit
        + area_exit
        * (pressure_exit - pressure_ambient)
        / mass_flowrate
    )
    thrust = constants.g * mass_flowrate * specific_impulse

def velocity_sound(t: float, rdot: float, y: float):
    """Velocity of Sound"""
    return np.sqrt(y * rdot * t)

def velocity_from_mach(temperature, mach, ratio_specific_heat, mass_molar):
    """Velocity from Mach"""
    return mach * np.sqrt(
        ratio_specific_heat * constants.R * temperature / mass_molar
    )

def density_from_mach(pressure, temperature, mass_molar):
    """Density from Mach"""
    return pressure / ((constants.R / mass_molar) * temperature)

def temp_from_mach(mach, temperature_chamber, ratio_specific_heat):
    """Temp from Mach"""
    return temperature_chamber / (
        1 + (ratio_specific_heat - 1) / 2 * mach**2
    )

def pressure_from_mach(mach, pressure_chamber, ratio_specific_heat):
    """Pressure from mach"""
    return pressure_chamber / (
        1 + (ratio_specific_heat - 1) / 2 * mach**2
    ) ** (ratio_specific_heat / (ratio_specific_heat - 1))

def mass_flux(mach, ratio_specific_heat, pressure_chamber, temperature_chamber, mass_molar):
    """Mass Flux from Mach"""
    temperature = temp_from_mach(mach, temperature_chamber, ratio_specific_heat)
    k = ratio_specific_heat
    return (
        pressure_chamber
        / temperature
        * (temperature / temperature_chamber) ** (k / (k - 1))
        * np.sqrt(
            2
            * mass_molar
            * (temperature_chamber - temperature)
            / (constants.R * (k - 1))
        )
    )

def calc_area_ratio(mach, ratio_specific_heat):
    """Area Ratio"""
    k = ratio_specific_heat
    return 1 / mach * ((2 + (k - 1) * mach**2) / (k + 1)) ** ((k + 1) / (2 * k - 2))

def area_ratio_derivative(mach, ratio_specific_heat):
    """Area Ratio Derivative"""
    k = ratio_specific_heat
    return (
        (2 * mach**2 - 2)
        * (((k - 1) * mach**2 + 2) / (k + 1)) ** ((k + 1) / (2 * k - 2))
        / ((k - 1) * mach**4 + 2 * mach**2)
    )

def mach_from_area(area, area_throat, ratio_specific_heat, is_subsonic: bool = False):
    """Mach from area"""

    max_iterations = 32
    threshold = 10**-16
    area_ratio = area / area_throat

    if is_subsonic:
        mach_guess = 1 / area_ratio
        for i in range(max_iterations):
            error = calc_area_ratio(mach_guess, ratio_specific_heat) - area_ratio
            if all(abs(error)) >= threshold:
                slope = area_ratio_derivative(mach_guess, ratio_specific_heat)
                mach_guess -= error / slope
    else:
        a = 2 / (ratio_specific_heat + 1)
        mach_guess = 1 + np.sqrt(abs(area_ratio - 1)) / np.sqrt(a)
        for i in range(max_iterations):
            error = calc_area_ratio(mach_guess, ratio_specific_heat) - area_ratio
            if all(abs(error)) >= threshold:
                slope = area_ratio_derivative(mach_guess, ratio_specific_heat)
                mach_guess -= error / slope

    return mach_guess