import matplotlib as mpl

# 5 Feb 2022
# https://stackoverflow.com/questions/60870584/importing-siunitx-into-python
#mpl.use("pgf")

# 5 Feb 2022
# https://matplotlib.org/stable/gallery/userdemo/pgf_preamble_sgskip.html
# https://tex.stackexchange.com/questions/391074/how-to-use-the-siunitx-package-within-python-matplotlib
# https://stackoverflow.com/questions/11354149/python-unable-to-render-tex-in-matplotlib
# https://matplotlib.org/stable/tutorials/text/pgf.html
mpl.rcParams.update({
    "figure.figsize": (9, 6),
    "figure.dpi": 200,
    "font.family": "serif",
    "font.size": 20,
    "text.usetex": True,
    'text.latex.preamble': "\n".join([
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{bm}",
#        r"\usepackage{siunitx}",
        ])
    #"pgf.rcfonts": False,
    #"pgf.texsystem": "pdflatex",
    #"pgf.preamble": "\n".join([
        #r"\usepackage{amsmath}",
        #r"\usepackage[utf8]{inputenc}",
        #r"\usepackage[T1]{fontenc}",
        #r"\usepackage{siunitx}", # load additional packages
    #])
})


legend_dict = {"bbox_to_anchor": (0.5, 1.3), "ncol": 4, "loc": "upper center"}
