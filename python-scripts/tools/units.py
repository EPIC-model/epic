units = {"time": "s", "position": "m", "velocity": "m/s", "energy": "$m^4/s^2$"}


def get_label(label, unit):
    if unit == None:
        return label
    else:
        return label + " (" + unit + ")"
