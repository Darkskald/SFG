import numpy as np


# todo: standardize verbs in function names like "calc", "yield", "get"

def calculate_area(diameter: float):
    return np.pi * (diameter / 2) ** 2


def calculate_molecule_count(volume: float, concentration: float):
    volume_in_liter = volume / 1e6
    amount_in_milli_mole = volume_in_liter * concentration
    total_amount = 6e23 * 1e-3 * amount_in_milli_mole
    return total_amount


def calculate_area_per_molecule(volume, concentration, diameter=5.2):
    amount = calculate_molecule_count(volume, concentration)
    area = calculate_area(diameter) * 1e16
    return area / amount
