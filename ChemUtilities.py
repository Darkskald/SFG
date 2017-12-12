#calculate amount of substance


def amount_from_mass(mass, molar_mass):
    amount = (mass/molar_mass)
    return amount


def amount_from_conc(conc, volume):
    amount = conc* volume
    return amount


#c1*c2=c2*v2

def calc_v1(c2, v2, c1):
    """Dilution calculation"""
    v1 = (c2*v2)/c1
    return v1


def calc_c2(c1, v1, v2):
    """Titration"""
    c2 = (c1*v1)/v2
    return c2

