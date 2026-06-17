try:
    from scipy.constants import e
except ImportError:
    e = 1.60217663e-19



TIME = {
    "s":  1.0,
    "ms": 1.e-3,
    "µs": 1.e-6,
    "ns": 1.e-9
    }

LENGTH = {
    "km": 1.e3,
    "m":  1.0,
    "cm": 1.e-2,
    "mm": 1.e-3,
    "µm": 1.e-6,
    "nm": 1.e-9,
    "pm": 1.e-12,
    "fm": 1.e-15
}

MASS = {
    "g":  1.0,
    "mg": 1.e-3,
    "kg": 1.e3
    }

ENERGY = {
    "J":  1.0,
    "Ws": 1.0,
    "eV": e
    }



QUANTITIES = {
    "time":   TIME,
    "length": LENGTH,
    "mass":   MASS,
    "energy": ENERGY
    }



# SI base unit equivalents
SI = {
    "force":         ("N",  "kg * m * s**(-2)"),
    "pressure":      ("Pa", "kg * m**(-1) * s**(-2)"),
    "energy":        ("J",  "kg * m**2 * s**(-2)"),
    "power":         ("W",  "kg * m**2 * s**(-3)"),
    "magnetic flux": ("Wb", "kg * m**2 * s**(-2) * A**(-1)")
    }
