units = {
    'time'      : 's',
    'position'  : 'm',
    'velocity'  : 'm/s'
}


def get_label(label, unit):
    if unit == None:
        return label
    else:
        return label + ' (' + unit + ')'
