#!/usr/bin/env python
# ======================================================================
#              Merges two ellipses in two different ways
# ======================================================================

# =====perform various generic imports=====
import warnings, os, sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc

rcParams.update({"figure.autolayout": True})
warnings.simplefilter("ignore", DeprecationWarning)

## global settings

# set tick label size:
label_size = 25
mpl.rcParams["xtick.labelsize"] = label_size
mpl.rcParams["ytick.labelsize"] = label_size
# set x tick width and size:
mpl.rcParams["xtick.major.size"] = 10
mpl.rcParams["xtick.major.width"] = 2
mpl.rcParams["xtick.minor.size"] = 5
mpl.rcParams["xtick.minor.width"] = 1
# set y tick width and size:
mpl.rcParams["ytick.major.size"] = 10
mpl.rcParams["ytick.major.width"] = 2
mpl.rcParams["ytick.minor.size"] = 5
mpl.rcParams["ytick.minor.width"] = 1

# Ensure latex fonts throughout:
rc("font", **{"family": "Times New Roman"})
rc("text", usetex=True)
# =========================================

if "-h" in sys.argv or "--help" in sys.argv:
    print("Merges two ellipses in two different ways")
    exit(0)

# ---------------------------------------------------------------
# Input all parameters:
print()
print(" The first ellipsoid is centred at the origin in standard position.")
print(" We take its area to be pi.")
print()
x1 = 0.0
y1 = 0.0
a1b1 = 1.0
asp1_in = input(" Enter a_1/b_1 (default 2): ")
asp1 = float(asp1_in or 2.0)
# B_1 matrix:
b111 = asp1
b112 = 0.0
b122 = 1.0 / asp1

print()
print(" The second ellipsoid is centred at (X_2,Y_2), has axes a_2 & b_2,")
print(" and is rotated by an angle phi_2 w.r.t. the positive x axis.")
print()
a2b2_in = input(" Enter a_2*b_2 (default 0.25): ")
a2b2 = float(a2b2_in or 0.25)
asp2_in = input(" Enter the aspect ratio a_2/b_2 (default 3): ")
asp2 = float(asp2_in or 3.0)
phi2_in = input(" Enter phi_2 (degrees, default 45): ")
phi2 = float(phi2_in or 45.0)
phi2 = phi2 * np.pi / 180.0
cphi2 = np.cos(phi2)
sphi2 = np.sin(phi2)
c2phi2 = np.cos(2.0 * phi2)
s2phi2 = np.sin(2.0 * phi2)
# B_2 matrix:
fp = 0.5 * a2b2 * (asp2 + 1.0 / asp2)
fm = 0.5 * a2b2 * (asp2 - 1.0 / asp2)
b211 = fp + fm * c2phi2
b212 = fm * s2phi2
b222 = fp - fm * c2phi2
print()
print(" Let X_2 = R*cos(theta) and Y_2 = R*sin(theta).")
print(" Also, let r^2 = a*b = a_1*b_1 + a_2*b_2.")
print()
rrat_in = input(" Enter R/r (default 2): ")
rrat = float(rrat_in or 2.0)
theta_in = input(" Enter theta (degrees, default 30): ")
theta = float(theta_in or 30.0)
theta = theta * np.pi / 180.0

# ---------------------------------------------------------------
#
chi = 2.0 * rrat * np.cos(theta)
eta = 2.0 * rrat * np.sin(theta)

ab = a1b1 + a2b2
r = np.sqrt(ab)
# (X_2,Y_2):
x2 = 0.5 * r * chi
y2 = 0.5 * r * eta

# mu_1 & mu_2:
m1 = a1b1 / ab
m2 = a2b2 / ab

m1s = m1 * m1
m1m2 = m1 * m2
m2s = m2 * m2

# B* matrix:
bs11 = m1m2 * chi * chi + m1s * b111 + m2s * b211
bs12 = m1m2 * chi * eta + m1s * b112 + m2s * b212
bs22 = m1m2 * eta * eta + m1s * b122 + m2s * b222

# Solve a quartic to find best fit ellipse:
detbs = bs11 * bs22 - bs12 * bs12
coeffs = np.array(
    [1.0, 0.0, -2.0 - detbs, bs11 ** 2 + bs22 ** 2 + 2.0 * bs12 ** 2, 1.0 - detbs]
)
q = np.roots(coeffs)

mu = -coeffs[4] / coeffs[3]
merr = 1.0
while merr > 1.0e-12:
    mup = -(coeffs[4] + mu * (coeffs[3] + mu * (coeffs[2] + mu * mu))) / (
        coeffs[3] + mu * (2.0 * coeffs[2] + 4.0 * mu * mu)
    )
    mu = mu + mup
    merr = abs(mup)

print(" Newton-Raphsen iteration gives mu = %12.9f" % mu)

print()
print(" The real roots of the quartic and the F norms are as follows:")
print()
print("      mu              F")
fmin = 1.0e20
for k in range(4):
    if abs(np.imag(q[k])) < 1.0e-9:
        mu = np.real(q[k])
        b11 = (bs11 - mu * bs22) / (1.0 - mu ** 2)
        b22 = (bs22 - mu * bs11) / (1.0 - mu ** 2)
        b12 = bs12 / (1.0 - mu)
        f = 0.5 * (b11 - bs11) ** 2 + (b12 - bs12) ** 2 + 0.5 * (b22 - bs22) ** 2
        print(" %12.9f   %12.9f" % (mu, f))
        if f < fmin:
            fmin = f
            mumin = mu

print()
mu = mumin
print(" Selecting mu = %12.9f" % mu)

b11 = (bs11 - mu * bs22) / (1.0 - mu ** 2)
b22 = (bs22 - mu * bs11) / (1.0 - mu ** 2)
b12 = bs12 / (1.0 - mu)

# (X,Y):
x = m1 * x1 + m2 * x2
y = m1 * y1 + m2 * y2

# a/b:
bsum = 0.5 * (b11 + b22)
bdif = 0.5 * (b11 - b22)
asp = bsum + np.sqrt(bdif ** 2 + b12 ** 2)
print()
print(" Aspect ratio a/b of merged ellipse = %12.9f" % asp)
print()

# ahat vector or (cos(phi),sin(phi)):
adif = asp - b22
den = np.sqrt(adif ** 2 + b12 ** 2)
cphi = adif / den
sphi = b12 / den

# Next work out characteristics of "geometric" ellipse found by scaling B*:
bmult = 1.0 / np.sqrt(detbs)
bg11 = bs11 * bmult
bg12 = bs12 * bmult
bg22 = bs22 * bmult
bsum = 0.5 * (bg11 + bg22)
bdif = 0.5 * (bg11 - bg22)
aspg = bsum + np.sqrt(bdif ** 2 + bg12 ** 2)
print(" Aspect ratio a/b of geometric ellipse = %12.9f" % aspg)
print()
adif = aspg - bg22
den = np.sqrt(adif ** 2 + bg12 ** 2)
cphig = adif / den
sphig = bg12 / den

f = 0.5 * (bg11 - bs11) ** 2 + (bg12 - bs12) ** 2 + 0.5 * (bg22 - bs22) ** 2
print(" Error in fit = %12.9f" % f)

# ------------------------------------------------------------------------
# Set up figure:
fig1 = plt.figure(1, figsize=[6, 6])
ax1 = fig1.add_subplot(111)
ax1.set_xlabel("$x$", fontsize=30)
ax1.set_ylabel("$y$", fontsize=30)

alpha = np.linspace(0.0, 2.0 * np.pi, 361)
calp = np.cos(alpha)
salp = np.sin(alpha)

# Plot original ellipses in red:
a1 = np.sqrt(a1b1 * asp1)
b1 = np.sqrt(a1b1 / asp1)
xp = x1 + a1 * calp
yp = y1 + b1 * salp
ax1.plot(xp, yp, "r", lw=2)
xy1 = max(np.amax(abs(xp)), np.amax(abs(yp)))

a2 = np.sqrt(a2b2 * asp2)
b2 = np.sqrt(a2b2 / asp2)
xt = a2 * calp
yt = b2 * salp
xp = x2 + xt * cphi2 - yt * sphi2
yp = y2 + yt * cphi2 + xt * sphi2
ax1.plot(xp, yp, "r", lw=2)
xy2 = max(np.amax(abs(xp)), np.amax(abs(yp)))

# Plot merged ellipse in black:
a = np.sqrt(ab * asp)
b = np.sqrt(ab / asp)
xt = a * calp
yt = b * salp
xp = x + xt * cphi - yt * sphi
yp = y + yt * cphi + xt * sphi
ax1.plot(xp, yp, "k", lw=2, label="optimal")
xy = max(np.amax(abs(xp)), np.amax(abs(yp)))

# Plot geometric ellipse in magenta (dashed):
ag = np.sqrt(ab * aspg)
bg = np.sqrt(ab / aspg)
xt = ag * calp
yt = bg * salp
xp = x + xt * cphig - yt * sphig
yp = y + yt * cphig + xt * sphig
ax1.plot(xp, yp, "m--", lw=2, label="geometric")
xyg = max(np.amax(abs(xp)), np.amax(abs(yp)))

ax1.legend(loc="lower left", prop={"size": 20})

# Determine a nice plot window:
xymax = 1.05 * max(xy1, xy2, xy, xyg)
ax1.set_xlim(-xymax, xymax)
ax1.set_ylim(-xymax, xymax)

plt.show()
