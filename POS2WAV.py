#Mapping the polar coordinates of a static wave source to the phase differences upon initial reception.

import numpy as np
class Receivers:
    def __init__(self, dimension, shape, angles, sides):
        self.form = shape
        self.dim = dimension
        if self.form == "triangle":
            self.ang_k = angles[0] # angle of vertex at (0,0)
            self.sdlen = sides[:2]

class Wave_func:
    def __init__(self, frequency, wavelength, amplitude):
        self.f = frequency
        self.lam = wavelength
        self.amp = amplitude

class Initial_wave_info:
    def __init__(self, frequency, wavelength, amplitude, distance_differences, time_differences, phase_differences):
        self.f = frequency
        self.lmbd = wavelength
        self.spd = frequency * wavelength
        self.amp = amplitude
        self.dd = distance_differences
        self.td = time_differences
        self.phd = phase_differences
    def __repr__(self):
        return 'Source to verteces \ndistance differences: p0,1 {0}; p0,2 {1} \ntime differences: p0,1 {2}; p0,2 {3} \nphase differences: p0,1 {4}; p0,2 {5}'.format(self.dd[0], self.dd[1], self.td[0], self.td[1], self.phd[0], self.phd[1])

test_tri_recv = Receivers(2, "triangle", [np.pi/2], [1, 1])

# map wave source position to receivers' phase differences, 2D
# polar coordinates

def tri_map_p2w(recv, dist_0, theta_0, freq, wvlen, ampl):
    r_01 = recv.sdlen[0]
    r_02 = recv.sdlen[1]
    theta_1 = np.pi - theta_0
    theta_2 = 2*np.pi - theta_1 - recv.ang_k
    d_diff = []
    t_diff = []
    ph_diff = []
    if theta_1 == np.pi/2:
        x_01 = 0 #distance difference between p0 (0,0) and p1 (r_01, pi)
    if theta_1 < np.pi/2: 
        x_01 = ((-2)*dist_0 - np.sqrt(4*dist_0**2 - 4*(2*np.cos(theta_1)*dist_0*r_01 - r_01**2)))/2 #consult documentation
    if theta_1 > np.pi/2: 
        x_01 = ((-2)*dist_0 + np.sqrt(4*dist_0**2 - 4*(2*np.cos(theta_1)*dist_0*r_01 - r_01**2)))/2
    d_diff.append(x_01)
    if theta_2 == np.pi/2:
        x_02 = 0 
    if theta_2 < np.pi/2: 
        x_02 = ((-2)*dist_0 - np.sqrt(4*dist_0**2 - 4*(2*np.cos(theta_2)*dist_0*r_02 - r_02**2)))/2 #consult documentation
    if theta_2 > np.pi/2: 
        x_02 = ((-2)*dist_0 + np.sqrt(4*dist_0**2 - 4*(2*np.cos(theta_2)*dist_0*r_02 - r_02**2)))/2
    d_diff.append(x_02)
    t_diff_01 = x_01 / (freq * wvlen)
    t_diff.append(t_diff_01)
    t_diff_02 = x_02 / (freq * wvlen)
    t_diff.append(t_diff_02)
    ph_diff_01 = t_diff_01 * freq
    ph_diff.append(ph_diff_01)
    ph_diff_02 = t_diff_02 * freq
    ph_diff.append(ph_diff_02)
    init_w_info = Initial_wave_info(freq, wvlen, ampl, d_diff, t_diff, ph_diff)
    return init_w_info

test_wave = Wave_func(1000, 0.0005, 1)
test_map = tri_map_p2w(test_tri_recv, 200, np.pi/3, test_wave.f, test_wave.lam, test_wave.amp)
print(test_map)