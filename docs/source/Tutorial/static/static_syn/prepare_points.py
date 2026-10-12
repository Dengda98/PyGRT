import numpy as np

theta = np.linspace(0.0, 2.0 * np.pi, 72, endpoint=False)
radius = 5.0
points = np.column_stack((
    radius * np.cos(theta),
    radius * np.sin(theta),
    np.zeros_like(theta),
))
np.savetxt("rcv_pts.txt", points, fmt="%.8f", header="north east depth (km)")
