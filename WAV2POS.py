import numpy as np
import POS2WAV as p2w

def tri_map_w2p(recv, init_rec_w_info):
    ph_diff = init_rec_w_info.phd
    t_diff = init_rec_w_info.td
    d_diff = init_rec_w_info.dd
    #Key is to solve for dist_0 and theta_0 (using theta_1) with the equation set:
    #(d_0 + x_01)^2 = d_1^2 + r_01^2 - 2*r_01*d_1*cos(theta_1)      This comes from the law of cosines; consult documentation for detailed solving steps
    #(d_0 + x_02)^2 = d_1^2 + r_02^2 - 2*r_02*d_1*cos(theta_2)      (theta_2 = 2*pi - ang_k[0] - theta_1
    r_01 = recv.sdlen[0]
    r_02 = recv.sdlen[1]
    x_01 = d_diff[0]
    x_02 = d_diff[1]
    theta_k = recv.ang_k
    #Simplfying Params, consult documentation
    A_ = 2 * r_01 * (r_02**2 - x_02**2)
    B_ = 2 * r_02 * (r_01**2 - x_01**2)
    C_ = 2 * x_01 * (r_02**2 - x_02**2)
    D_ = 2 * x_02 * (r_01**2 - x_01**2)
    E_ = np.cos(2*np.pi - theta_k)
    F_ = np.sin(2*np.pi - theta_k)
    G_ = -2 * ((D_ - C_) * A_ - (D_ - C_) * B_ * E_ + A_ * B_ * E_)
    H_1 = (A_**2 + (B_ * E_)**2 - (B_ * F_)**2) #theta_1 range 0 to pi / theta_0 is positive / sin(theta_1) > 0
    H_2 = (A_**2 + (B_ * E_)**2 + (B_ * F_)**2) #theta_1 range pi to 2*pi / theta_0 is negative / sin(theta_1) < 0
    I_1 = ((D_ - C_)**2 + (B_ * F_)**2)
    I_2 = ((D_ - C_)**2 - (B_ * F_)**2)
    J_1 = G_**2 - 4 * H_1 * I_1
    J_2 = G_**2 - 4 * H_2 * I_2

    list_of_possible_cos_t1s = []
    if J_1 >= 0: #solution has to be real
        cos_theta_1_gr_pos = (-G_ + np.sqrt(J_1)) / (2 * H_1)
        list_of_possible_cos_t1s.append([cos_theta_1_gr_pos, 'pos'])
        cos_theta_1_sm_pos = (-G_ - np.sqrt(J_1)) / (2 * H_1)
        list_of_possible_cos_t1s.append([cos_theta_1_sm_pos, 'pos'])
    if J_2 >= 0:
        cos_theta_1_gr_neg = (-G_ + np.sqrt(J_2)) / (2 * H_2)
        list_of_possible_cos_t1s.append([cos_theta_1_gr_neg, 'neg'])
        cos_theta_1_sm_neg = (-G_ - np.sqrt(J_2)) / (2 * H_2)
        list_of_possible_cos_t1s.append([cos_theta_1_sm_neg, 'neg'])

    list_of_possible_t1s = []
    for ct1 in list_of_possible_cos_t1s:
        if ct1[1] == 'pos': #whether the result from arccos has to be adjusted
            list_of_possible_t1s.append(np.arccos(ct1[0]))
        else:
            list_of_possible_t1s.append(2 * np.pi - np.arccos(ct1[0]))
    list_of_possible_t0s = []
    list_of_possible_d0s = []
    ret_list = []
    for t1 in list_of_possible_t1s:
        list_of_possible_t0s.append(np.pi - t1)
        list_of_possible_d0s = (r_01**2 - x_01**2) / (2 * x_01 + 2 * r_01 * np.cos(t1))
        ret_list.append('possible theta_0: {0}, correspondent dist_0: {1}'.format(np.pi - t1, (r_01**2 - x_01**2) / (2 * x_01 + 2 * r_01 * np.cos(t1))))

    return ret_list
    #To select the correct case among 4, consider the signs of x_01 and x_02
    #The wave source point falls on one of the sections of space divided by the lines extended from the perpendicular line segment bisectors of each side of the triangle

print(tri_map_w2p(p2w.test_tri_recv, p2w.test_map))