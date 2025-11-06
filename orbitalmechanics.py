import requests
from sgp4.api import Satrec, jday
from datetime import datetime, timedelta
import numpy as np
from scipy.spatial.distance import cdist
from poliastro.plotting.static import StaticOrbitPlotter
from poliastro.twobody.orbit import Orbit
from poliastro.bodies import Earth
from astropy import units as u
import matplotlib.pyplot as plt

# TLE for Satellite Starlink-11106
line1 = "1 59765C 24091N   25306.82340278 -.00007572  00000+0 -62396-4 0  3062"
line2 = "2 59765  53.1560 111.3788 0001142  82.3215 200.8315 15.69709430    16"

sat = Satrec.twoline2rv(line1, line2)

# Position computation 
jd, fr = jday(2025, 11, 3, 12, 0, 0)
e, r, v = sat.sgp4(jd, fr)

print("Position (km):", r)
print("Velocity (km/s):", v)