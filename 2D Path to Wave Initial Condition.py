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
        return 'Source to verteces \ndistance differences: p1,2 {0}; p1,3 {1} \ntime differences: p1,2 {2}; p1,3 {3} \nphase differences: p1,2 {4}; p1,3 {5}'.format(self.dd[0], self.dd[1], self.td[0], self.td[1], self.phd[0], self.phd[1])

test_tri_recv = Receivers(2, "triangle", [np.pi/2], [1, 1])

# map wave source position to receivers' phase differences, 2D
# polar coordinates

def tri_map_p2w(recv, dist_1, theta_1, freq, wvlen, ampl):
    r_12 = recv.sdlen[0]
    r_13 = recv.sdlen[1]
    theta_2 = np.pi - theta_1
    theta_3 = 2*np.pi - theta_2 - recv.ang_k
    d_diff = []
    t_diff = []
    ph_diff = []
    if theta_2 == np.pi/2:
        x_12 = 0 #distance difference between p1 (0,0) and p2 (r_12, pi)
    if theta_2 < np.pi/2: 
        x_12 = ((-2)*dist_1 - np.sqrt(4*dist_1**2 - 4*(2*np.cos(theta_2)*dist_1*r_12 - r_12**2)))/2 #consult documentation
    if theta_2 > np.pi/2: 
        x_12 = ((-2)*dist_1 + np.sqrt(4*dist_1**2 - 4*(2*np.cos(theta_2)*dist_1*r_12 - r_12**2)))/2
    d_diff.append(x_12)
    if theta_3 == np.pi/2:
        x_13 = 0 
    if theta_3 < np.pi/2: 
        x_13 = ((-2)*dist_1 - np.sqrt(4*dist_1**2 - 4*(2*np.cos(theta_3)*dist_1*r_13 - r_13**2)))/2 #consult documentation
    if theta_3 > np.pi/2: 
        x_13 = ((-2)*dist_1 + np.sqrt(4*dist_1**2 - 4*(2*np.cos(theta_3)*dist_1*r_13 - r_13**2)))/2
    d_diff.append(x_13)
    t_diff_12 = x_12 / (freq * wvlen)
    t_diff.append(t_diff_12)
    t_diff_13 = x_13 / (freq * wvlen)
    t_diff.append(t_diff_13)
    ph_diff_12 = t_diff_12 * freq
    ph_diff.append(ph_diff_12)
    ph_diff_13 = t_diff_13 * freq
    ph_diff.append(ph_diff_13)
    init_w_info = Initial_wave_info(freq, wvlen, ampl, d_diff, t_diff, ph_diff)
    return init_w_info

test_wave = Wave_func(1000, 0.0005, 1)
test_map = tri_map_p2w(test_tri_recv, 200, np.pi/3, test_wave.f, test_wave.lam, test_wave.amp)
print(test_map)