from dynamics import *
import pickle
import os

bd_step_size_fs = 1000.0  # simulation time step in femotoseconds (10^-15 sec)


def create_pickle_from_array(np_array, file_name, save_path):
    """
    method to save pickles array
    """
    pickle.dump(np_array,
                open(os.path.join(save_path, file_name + ".pkl"), "wb"))


def create_model(seq, nchains, rmf_filename, bead_radius, sphere_radius, kbs,
                 nres_per_bead, k_in, k_out, is_center=True):
    """## Building a dynamic model - parts, interactions, dynamics
    OK, let's begin by building our first dynamic model!
    As we learned in class, a dynamic model stands on three pillars:
    parts, interactions, and dynamics.

    But first we need to construct an [IMP `model`]
    (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1Model.html):
    """

    m = IMP.Model()

    """### I. Add parts 
    #### Generate strings of beads
    In IMP, model parts are called [Particles]
    (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1Particle.html).
    The beads in our `ProteinChain` object are such particles - in this case, 
    they represent physical objects with Cartesian coordinates. 

    To generate two string of beads, we use the `ProteinChainFactory` class
    `create()` method:  
    """

    protein_chain_factory = ProteinChainFactory \
        (model=m,  # the model in which the beads reside
         default_radius_A=bead_radius,  # radius of a bead
         k_in=k_in,
         k_out=k_out,
         # the force coefficient for the spring holding consecutive beads together
         # (large number = stiff spring)
         relative_rest_distance=3.0,
         # the resting distance between bead centers, relative to the radius of single bead
         nres_per_bead=nres_per_bead,
         kbs=kbs,

         sphere_radius=sphere_radius)  # number of residues per bead in the string-of-beads

    label = "popZ"
    # nchains = 2
    chains = []
    for i in range(nchains):
        chain = protein_chain_factory.create(seq, f"{label}_{i}", is_center)
        chains.append(chain)

    """### Keep chains in a hierarchical data structure

    Next, we create a hierarchical data structure, in which the two chain descend
     from a common root particle called `p_root`. Note that ProteinChainFactory made sure that the beads of each chain also descend from a common root `protein_chain.root_as_h`.



    """

    p_root = IMP.Particle(m, "root")
    h_root = IMP.atom.Hierarchy.setup_particle(p_root)  # decorator
    for chain in chains:
        h_root.add_child(chain.root_as_h)

    """**Explanation**: IMP Particles don't have to represent actual physical 
    particles. To allow *particle* `p_root` to have children,
    we [decorate](https://https://en.wikipedia.org/wiki/Decorator_pattern)
    it using the [`IMP.atom.Hierarchy` class]
    (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1atom_1_1Hierarchy.html).
    The [decorator design pattern](https://https://en.wikipedia.org/wiki/Decorator_pattern)
     adds new functions to existing classes at runtime instead of during compile time. Think about `h_root` as particle `p_root` cast as class Hierarchy.

    ### II. Add interactions
    Interactions in IMP are created using the concept of a "restraint" acting on a group of particles (a single particle, a pair of particles, etc.). All restraints derive from the class [`Restraint`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1Restraint.html). 

    We first define an [`ExcludedVolumeRestraint`]
    (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1core_1_1ExcludedVolumeRestraint.html)
    object - telling IMP that different beads should not overlap with each other, 
    at least not too much:
    """

    # Add excluded volume restraints among all (close pairs of) particles:
    evr = IMP.core.ExcludedVolumeRestraint(get_all_beads(chains),
                                           # particles to be restraints
                                           1.0,  # force constant
                                           10.0,
                                           # slack parameter affecting speed only
                                           "Excluded-Volume"
                                           # a string identifier
                                           )

    """The springs between the beads of each protein chain were already defined by 
    ProteinChainFactory, so now we create `rsf`, a scoring function 
    (a.k.a. energy function) made of our restraints."""

    restraints = [chain.restraint for chain in chains] + [evr]
    for chain0 in chains:
        for chain1 in chains:
            if chain0 != chain1:
                restraints.append(
                    protein_chain_factory.get_interchain_restraint(chain0,
                                                                   chain1))
    restraints.append(protein_chain_factory.get_bounding_sphere_restraint(
        get_all_beads(chains)))

    rsf = IMP.core.RestraintsScoringFunction(restraints,
                                             "Scoring function")  # Energy function

    """### III. Add dynamics
    Finally, we create a [`BrownianDynamics`]
    (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1atom_1_1BrownianDynamics.html)
     object, which enables us to simulate our interacting parts:
    """
    bd = start_bd_simulation(rsf, m)
    T_ns, E, D, chains_on_iteration = create_trajectory_file(rmf_filename, bd,
                                                             h_root,
                                                             restraints, m,
                                                             chains,
                                                             sphere_radius,
                                                             k_out)
    return T_ns, E, D, chains_on_iteration


def start_bd_simulation(rsf, model):
    # BD
    bd = IMP.atom.BrownianDynamics(model)
    bd.set_scoring_function(rsf)  # tell BD about our scoring (energy) function
    bd.set_maximum_time_step(bd_step_size_fs)
    T = 300  # temperature in Kalvin
    bd.set_temperature(T)
    return bd


def create_trajectory_file(rmf_filename, bd, h_root, restraints, model, chains,
                           sphere_radius, k_out):
    """#### Making a movie (trajectory) file
    We conclude by telling our simulation to output a trajectory (movie)
    file every 10000 frames. The movie is saved using the [Rich Molecular File format]
    (https://integrativemodeling.org/rmf/format.html), which can be viewed using e.g. the [Chimera software](https://www.cgl.ucsf.edu/chimera/download.html) or [ChimeraX](https://cxtoolshed.rbvi.ucsf.edu/apps/chimeraxrmf).

    **Note:** if you get a "UsageError: "Opening a file that is still being written
    is asking for trouble." error for any reson, just change the filename variable
    """

    rmf = RMF.create_rmf_file(rmf_filename)
    rmf.set_description("Brownian dynamics trajectory with {}fs timestep.\n" \
                        .format(bd_step_size_fs))
    IMP.rmf.add_hierarchy(rmf,
                          h_root)  # Telling the movie that it should save all descendents of h_root
    IMP.rmf.add_restraints(rmf,
                           restraints)  # the restraints are also saved and can be viewed in Chimera

    IMP.rmf.add_geometry(rmf, IMP.display.SphereGeometry(
        IMP.algebra.Sphere3D(IMP.algebra.Vector3D(0, 0, 0), sphere_radius)))

    sos = IMP.rmf.SaveOptimizerState(model,
                                     rmf)  # an optimizer state is invoked every n frames of simulation by the BD simulation
    sos.set_simulator(bd)
    frames_interval = 10000
    sos.set_period(frames_interval)
    bd.add_optimizer_state(sos)
    sos.update_always(
        "initial conformation")  # save the initial conformation to the RMF file

    """## Simulating dynamics
    Without further ado, let's simulate our system using `bd.optimize(n)` where `n` is the number of time steps to simulate. 

    We shall keep a few statistics every `n_inner_cycle` frames - the simulation time, the energy of the system (our scoring function), and the distance between the first and last bead of each chain.

    """

    T_ns = []  # time in nanoseconds
    E = []  # energy
    D = [[] for chain in chains]
    chains_on_iteration = []
    n_outer = 50000  # outer loop number of iterations
    n_inner = 250  # optimization per iteration

    for i in range(n_outer):
        time_fs = bd.get_current_time()
        time_ns = time_fs * 1e-6  # a nanosecond is a million femtoseconds
        if i % (n_outer // 10) == 0:
            print(f"Simulated for {time_ns:.1f} nanoseconds so far")
        bd.optimize(n_inner)
        T_ns.append(time_ns)  # keep time
        E.append(bd.get_last_score())  # keep energy
        for i, chain in enumerate(chains):  # keep distances for each chain
            distance = IMP.core.get_distance(IMP.core.XYZ(chain.beads[0]),
                                             IMP.core.XYZ(chain.beads[-1]))
            D[i].append(distance)
        chains_on_iteration.append(generate_vecs_of_beads(chains))

    print(f"FINISHED. Simulated for {time_ns:.1f} nanoseconds in total.")

    rmf.flush()  # make sure RMF file is properly saved
    T_ns = np.array(T_ns)  # this will be convenient later
    E = np.array(E)  # likewise
    D = np.array(D)
    C = np.array(chains_on_iteration)

    # pickling array

    cwd = os.getcwd()
    kout = f'{k_out:.3f}'
    spr_rad = f'{sphere_radius:.3f}'

    create_pickle_from_array(T_ns, "T_ns_" + "kout_" +
                             kout + "_" + "spr_rad_" + spr_rad, cwd)
    create_pickle_from_array(E, "Energy_" + "kout_" +
                             kout + "_" + "spr_rad_" + spr_rad, cwd)
    create_pickle_from_array(D, "Distance-e-to-e_" + "kout_" +
                             kout + "_" + "spr_rad_" + spr_rad, cwd)
    create_pickle_from_array(C, "chain_on_iterations_" + "kout_" +
                             kout + "_" + "spr_rad_" + spr_rad, cwd)

    return T_ns, E, D, chains_on_iteration


def get_all_beads(chains):
    beads_set = set().union(*[chain.beads for chain in chains])
    return list(beads_set)


def generate_vecs_of_beads(chains):
    """
    method to create XYZs of chain beads array
    C[i][j] = XYZ of chain #i bead #j
    """
    vecs_of_beads = []
    vecs_of_chains = []
    for chain in chains:
        for bead in chain.beads:
            vecs_of_beads.append(np.array(IMP.core.XYZ(
                bead).get_coordinates()))
        vecs_of_chains.append(vecs_of_beads)
    return vecs_of_chains

# bead_radius = 10.0
# sphere_radius = 50 * bead_radius
# IS_DEBUG_K = False
# kbs = 0.1
# k_in = 5.0
# k_out = 0.05
# if IS_DEBUG_K:
#     kbs = -kbs
# nchains = 2
# nres_per_bead = 20
# is_center = False
# rmf_filename = f"my_trajectory_br{bead_radius}_sr{sphere_radius}_kbs{kbs}_kin{k_in}_kout{k_out}_nchains{nchains}_nres_per_bead{nres_per_bead}.rmf"
# seq = "MSDQSQEPTMEEILASIRRIISEDDAPAEPAAEAAPPPPPEPEPEPVSFDDEVLELTDPI" \
#       "APEPELPPLETVGDIDVYSPPEPESEPAYTPPPAAPVFDRDEVAEQLVGVSAASAAASAF" \
#       "GSLSSALLMPKDGRTLEDVVRELLRPLLKEWLDQNLPRIVETKVEEEVQRISRGRGA"
# # seq = "MSDQSQEPTMEEILASIRRI"
#
# # T_ns, E, D, chains_on_iteration = create_model(seq, nchains, rmf_filename, bead_radius, sphere_radius, kbs, nres_per_bead, k_in, k_out, is_center)
# # x=0
