units = {"time": "s", "position": "m", "velocity": "m/s", "energy": "m^4/s^2"}


def get_label(label, unit, is_bokeh=False):
    if unit == None:
        return label
    else:
        if is_bokeh:
            return label + " $$(" + unit + ")$$"
        else:
            return label + " $(" + unit + ")$"
