import numpy as np
import constants

def performance(self):
    """
    Calculates performance parameters along the contour
    """

    self.contour_temperature = self.temp_from_mach(self.contour_mach)
    self.contour_pressure = self.pressure_from_mach(self.contour_mach)
    self.contour_density = self.density_from_mach(
        self.contour_pressure, self.contour_temperature
    )
    self.contour_mdot = self.mass_flux(self.contour_mach)
    self.contour_velocity = self.velocity_from_mach(
        temperature=self.contour_temperature, mach=self.contour_mach
    )
    self.specific_impulse = (
        self.velocity_exit
        + self.area_exit
        * (self.pressure_exit - self.pressure_ambient)
        / self.mass_flowrate
    )
    self.thrust = constants.g * self.mass_flowrate * self.specific_impulse

def velocity_sound(self, t: float, rdot: float, y: float):
    """Velocity of Sound"""
    return np.sqrt(y * rdot * t)

def velocity_from_mach(self, temperature, mach):
    """Velocity from Mach"""
    return mach * np.sqrt(
        self.ratio_specific_heat * constants.R * temperature / self.mass_molar
    )

def density_from_mach(self, pressure, temperature):
    """Density from Mach"""
    return pressure / ((constants.R / self.mass_molar) * temperature)

def temp_from_mach(self, mach):
    """Temp from Mach"""
    return self.temperature_chamber / (
        1 + (self.ratio_specific_heat - 1) / 2 * mach**2
    )

def pressure_from_mach(self, mach):
    """Pressure from mach"""
    return self.pressure_chamber / (
        1 + (self.ratio_specific_heat - 1) / 2 * mach**2
    ) ** (self.ratio_specific_heat / (self.ratio_specific_heat - 1))

def mass_flux(self, mach):
    """Mass Flux from Mach"""
    temperature = self.temp_from_mach(mach)
    k = self.ratio_specific_heat
    return (
        self.pressure_chamber
        / temperature
        * (temperature / self.temperature_chamber) ** (k / (k - 1))
        * np.sqrt(
            2
            * self.mass_molar
            * (self.temperature_chamber - temperature)
            / (constants.R * (k - 1))
        )
    )

def area_ratio(self, mach):
    """Area Ratio"""
    k = self.ratio_specific_heat
    return 1 / mach * ((2 + (k - 1) * mach**2) / (k + 1)) ** ((k + 1) / (2 * k - 2))

def area_ratio_derivative(self, mach):
    """Area Ratio Derivative"""
    k = self.ratio_specific_heat
    return (
        (2 * mach**2 - 2)
        * (((k - 1) * mach**2 + 2) / (k + 1)) ** ((k + 1) / (2 * k - 2))
        / ((k - 1) * mach**4 + 2 * mach**2)
    )

def mach_from_area(self, area, is_subsonic: bool = False):
    """Mach from area"""

    max_iterations = 32
    threshold = 10**-16
    area_ratio = area / self.area_throat

    if is_subsonic:
        mach_guess = 1 / area_ratio
        for i in range(max_iterations):
            error = self.area_ratio(mach_guess) - area_ratio
            if all(abs(error)) >= threshold:
                slope = self.area_ratio_derivative(mach_guess)
                mach_guess -= error / slope
    else:
        a = 2 / (self.ratio_specific_heat + 1)
        mach_guess = 1 + np.sqrt(abs(area_ratio - 1)) / np.sqrt(a)
        for i in range(max_iterations):
            error = self.area_ratio(mach_guess) - area_ratio
            if all(abs(error)) >= threshold:
                slope = self.area_ratio_derivative(mach_guess)
                mach_guess -= error / slope

    return mach_guess