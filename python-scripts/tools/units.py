units = {
    "buoyancy":     "m / s^2",
    "energy":       "m^4 / s^2",
    "position":     "m",
    "power":        "m^5 / s^4",
    "time":         "s",
    "velocity":     "m / s",
    "vorticity":    "1 / s",
    "wavenumber":   "1 / m"
}


def get_label(label, unit, is_bokeh=False):
    if unit == None:
        return label
    else:
        if is_bokeh:
            return label + r" $$(" + unit + ")$$"
        else:
            return label + r" (\si{" + unit + "})"
