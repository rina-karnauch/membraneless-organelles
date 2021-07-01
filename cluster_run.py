import numpy as np
import sys
import main_model

k_out_list = np.logspace(-3, -0.5, 15, base=10)
radius_list = np.logspace(2, 3, 8, base=10)

task_num = int(sys.argv[1])

sphere_radius = radius_list[task_num // 15]
k_out = k_out_list[task_num % 15]


# N_A = 6.022E23
# L_PER_A3 = 1E-27
#
# def get_sphere_vol(R):
#       return 4/3 * np.pi * (R ** 3)
#
# def get_concentration_M(n_chains, R):
#       return (n_chains / get_sphere_vol(R)) / N_A / L_PER_A3




# print(f"radius: {sphere_radius:.3f}, k_out: {k_out:.3f}")
#
# for k in k_out_list:
#       print(f"{k:.3f}")
# print("-------")
# for r in radius_list:
#       print(f"{r:.3f}")

# print(f"Concentrations in micromolar: {(1E6 * get_concentration_M(4, 100))}")

bead_radius = 10.0
# sphere_radius = radius
kbs = 0.1
k_in = 5.0
# k_out = 0.05

nchains = 2
nres_per_bead = 20
is_center = False
rmf_filename = f"my_trajectory_br{bead_radius}_sr{sphere_radius}_kbs{kbs}_kin{k_in}_kout{k_out}_nchains{nchains}_nres_per_bead{nres_per_bead}.rmf"
seq = "MSDQSQEPTMEEILASIRRIISEDDAPAEPAAEAAPPPPPEPEPEPVSFDDEVLELTDPI" \
      "APEPELPPLETVGDIDVYSPPEPESEPAYTPPPAAPVFDRDEVAEQLVGVSAASAAASAF" \
      "GSLSSALLMPKDGRTLEDVVRELLRPLLKEWLDQNLPRIVETKVEEEVQRISRGRGA"
# seq = "MSDQSQEPTMEEILASIRRI"

# T_ns, E, D, chains_on_iteration = create_model(seq, nchains, rmf_filename, bead_radius, sphere_radius, kbs, nres_per_bead, k_in, k_out, is_center)
# x=0
main_model.create_model(seq, nchains, rmf_filename, bead_radius, sphere_radius, kbs, nres_per_bead, k_in, k_out, is_center)

