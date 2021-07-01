import sys
import IMP
import IMP.atom
import IMP.algebra
import IMP.rmf
import IMP.core
import RMF
import IMP.container
import IMP.display
import IMP.npctransport
import numpy as np
import arviz as az
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm

sys.path.append('/usr/lib/python3.6/dist-packages')

# @title
from dataclasses import dataclass


@dataclass  # https://docs.python.org/3/library/dataclasses.html
class ProteinChain:
    """
    A string-of-beads representation of a protein chain
    from specified beads. The beads are connected by a restraint
    (think springs between consecutive beads)
    """
    # the "parent" of all beads in the chain
    root_p: IMP.Particle
    # the list of bead particles comprising the chain
    beads: list
    # a spring restraint on consecutive beads
    restraint: IMP.Restraint
    # The protein sequence that this chain represents
    sequence: str

    def __post_init__(
            self):  # this is called by automatically generated __init__
        # make sure root is the parent of all beads
        for bead in self.beads:
            assert IMP.atom.Hierarchy(bead).get_parent() == self.root_as_h

    @property
    def model(self) -> IMP.Model:
        return self.root_p.get_model()

    @property  # Hierarchy is a decorator
    def root_as_h(self) -> IMP.atom.Hierarchy:
        return IMP.atom.Hierarchy(self.root_p)


"""#### ProteinChainFactory class
We also included a class for generating `ProteinChain` objects for various
protein sequences. 
Details will be explained below, just run for now:
"""

# @title
FAKE_MASS = 1.0


class ProteinChainFactory:
    """
    A class for generating chains of beads
    """
    def __init__(self,
                 model: IMP.Model,
                 default_radius_A: float = 10.0,
                 # bead radius in angstroms (10^-10 m)
                 relative_rest_distance: float = 2.0,
                 # distance between bead centers, relative to radius
                 k_in: float = 1.0,
                 # force coefficient for spring connecting consecutive beads in kcal/mol/A^2 (force is k*distance
                 nres_per_bead: int = 5,
                 kbs: float = 0.1,
                 sphere_radius: float = 100,
                 k_out: float = 1.0
                 ):
        """
        :param model
        :param default_radius_A       Bead radius in angstroms
        :param relative_rest_distance Distance between bead centers, relative to radius
        :param k_kcal_per_mol_per_A2  Force coefficient for spring connecting consecutive beads
                                      in kcal/mol/A^2 units
        :param nres_per_bead          Number of protein residues represented by a single bead
        """
        self._model = model
        self._default_radius_A = default_radius_A
        self._rest_distance_A = relative_rest_distance * default_radius_A
        self._k_in = k_in
        self._k_out = k_out
        self._nres_per_bead = nres_per_bead
        self._kbs = kbs
        self._sphere_radius = sphere_radius
        self.s = IMP.algebra.Sphere3D(IMP.algebra.Vector3D(0, 0, 0), self._sphere_radius)
        self.init_locations = []

    @property
    def model(self): return self._model

    @property
    def default_radius_A(self): return self._default_radius_A

    @property
    def rest_distance_A(self): return self._rest_distance_A

    @property
    def k_in(self): return self._k_in

    @property
    def k_out(self): return self._k_out

    @property
    def nres_per_bead(self): return self._nres_per_bead

    def _create_bead(self,
                     name: str,
                     init_coord: IMP.algebra.Vector3D = None):
        p = IMP.Particle(self.model, name)
        p_as_xyzr = IMP.core.XYZR.setup_particle(
            p)  # A Decorator design pattern - adding functionality to an object at run time (~run-time inheritance)
        if init_coord is not None:
            p_as_xyzr.set_coordinates(init_coord)
        p_as_xyzr.set_coordinates_are_optimized(True)
        p_as_xyzr.set_radius(self.default_radius_A)
        IMP.atom.Mass.setup_particle(p, FAKE_MASS)  # required by Hierarchy
        IMP.atom.Diffusion.setup_particle(p)
        IMP.atom.Hierarchy.setup_particle(
            p)  # allow inclusion in IMP hierarchies
        IMP.display.Colored.setup_particle(p,
                                           IMP.display.get_display_color(0))
        return p

    def _create_restraint(self,
                          beads  # spring constant
                          ):
        hdps = IMP.core.HarmonicDistancePairScore(self.rest_distance_A,
                                                  self.k_in)
        cpc = IMP.container.ConsecutivePairContainer(self.model,
                                                     beads)  # convention - use abbreviation for multiword class names
        pr = IMP.container.PairsRestraint(hdps, cpc)
        return pr, hdps

    def get_interchain_restraint(self, chain0: ProteinChain,
                                 chain1: ProteinChain):
        center_A = 5
        threshold_A = 10
        thb = IMP.core.TruncatedHarmonicBound(center_A,
                                              self.k_out,
                                              threshold_A)
        dps = IMP.core.DistancePairScore(thb)
        pr = IMP.core.PairRestraint(chain0.model,
                                    dps,
                                    [chain0.beads[2], chain1.beads[2]])
        return pr

    def get_bounding_sphere_restraint(self, beads):
        f = IMP.core.HarmonicUpperBound(0, self._kbs) # mean = 0, k = 0.1 (k_kcal_per_mol_per_A2)
        bsss = IMP.core.BoundingSphere3DSingletonScore(f, self.s)
        r = IMP.container.SingletonsRestraint(bsss, beads)
        return r

    def create(self,
               sequence: str,  # protein sequence
               name: str,
               in_center: bool):  # a name of your choice for the protein chain
        p = IMP.Particle(self.model, name)
        p_as_h = IMP.atom.Hierarchy.setup_particle(
            p)  # allow inclusion in IMP hierarchies
        # add beads:
        n = len(sequence)
        nbeads = max(1, round(n / self.nres_per_bead))
        beads = []
        init_location = None
        if not in_center:
            init_location = self.get_free_location_in_sphere(self.init_locations)
        for i in range(nbeads):
            bead = self._create_bead(f"{name}_{i}", init_coord=init_location)
            p_as_h.add_child(bead)
            beads.append(bead)
        # restrain beads on a "string":
        restraint, harmonic_distance_pair_score = self._create_restraint(beads)

        return ProteinChain(root_p=p,
                            beads=beads,
                            restraint=restraint,
                            sequence=sequence)

    def get_free_location_in_sphere(self, ready_locations: [IMP.algebra.Vector3D]):
        new_location = IMP.algebra.get_random_vector_in(IMP.algebra.Sphere3D(IMP.algebra.Vector3D(0, 0, 0), self._sphere_radius/2))
        while not self.check_locations_collision(new_location, ready_locations, self.default_radius_A * 2):
            new_location = IMP.algebra.get_random_vector_in(self.s)
        return new_location

    def check_locations_collision(self,
                                  location_to_check: IMP.algebra.Vector3D,
                                  ready_locations: [IMP.algebra.Vector3D],
                                  threshold: float):
        for ready_location in ready_locations:
            if IMP.algebra.get_distance(location_to_check, ready_location) < threshold:
                return False
        return True

#
# """## Building a dynamic model - parts, interactions, dynamics
# OK, let's begin by building our first dynamic model!
# As we learned in class, a dynamic model stands on three pillars:
# parts, interactions, and dynamics.
#
# But first we need to construct an [IMP `model`]
# (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1Model.html):
# """
#
# m = IMP.Model()
#
# """### I. Add parts
# #### Generate strings of beads
# In IMP, model parts are called [Particles]
# (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1Particle.html).
# The beads in our `ProteinChain` object are such particles - in this case,
# they represent physical objects with Cartesian coordinates.
#
# To generate two string of beads, we use the `ProteinChainFactory` class
# `create()` method:
# """
#
# protein_chain_factory = ProteinChainFactory \
#     (model=m,  # the model in which the beads reside
#      default_radius_A=10.0,  # radius of a bead
#      k_kcal_per_mol_per_A2=5.0,
#      # the force coefficient for the spring holding consecutive beads together
#      # (large number = stiff spring)
#      relative_rest_distance=3.0,
#      # the resting distance between bead centers, relative to the radius of single bead
#      nres_per_bead=20)  # number of residues per bead in the string-of-beads
# seq = "MSDQSQEPTMEEILASIRRIISEDDAPAEPAAEAAPPPPPEPEPEPVSFDDEVLELTDPI" \
#       "APEPELPPLETVGDIDVYSPPEPESEPAYTPPPAAPVFDRDEVAEQLVGVSAASAAASAF" \
#       "GSLSSALLMPKDGRTLEDVVRELLRPLLKEWLDQNLPRIVETKVEEEVQRISRGRGA"
# label = "popZ"
# nchains = 2
# chains = []
# for i in range(nchains):
#     chain = protein_chain_factory.create(seq, f"{label}_{i}")
#     chains.append(chain)
#
# """### Keep chains in a hierarchical data structure
#
# Next, we create a hierarchical data structure, in which the two chain descend
#  from a common root particle called `p_root`. Note that ProteinChainFactory made sure that the beads of each chain also descend from a common root `protein_chain.root_as_h`.
#
#
#
# """
#
# p_root = IMP.Particle(m, "root")
# h_root = IMP.atom.Hierarchy.setup_particle(p_root)  # decorator
# for chain in chains:
#     h_root.add_child(chain.root_as_h)
#
# """**Explanation**: IMP Particles don't have to represent actual physical
# particles. To allow *particle* `p_root` to have children,
# we [decorate](https://https://en.wikipedia.org/wiki/Decorator_pattern)
# it using the [`IMP.atom.Hierarchy` class]
# (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1atom_1_1Hierarchy.html).
# The [decorator design pattern](https://https://en.wikipedia.org/wiki/Decorator_pattern)
#  adds new functions to existing classes at runtime instead of during compile time. Think about `h_root` as particle `p_root` cast as class Hierarchy.
#
# ### II. Add interactions
# Interactions in IMP are created using the concept of a "restraint" acting on a group of particles (a single particle, a pair of particles, etc.). All restraints derive from the class [`Restraint`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1Restraint.html).
#
# We first define an [`ExcludedVolumeRestraint`]
# (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1core_1_1ExcludedVolumeRestraint.html)
# object - telling IMP that different beads should not overlap with each other,
# at least not too much:
# """
#
#
# def get_all_beads(chains):
#     beads_set = set().union(*[chain.beads for chain in chains])
#     return list(beads_set)
#
#
# # Add excluded volume restraints among all (close pairs of) particles:
# evr = IMP.core.ExcludedVolumeRestraint(get_all_beads(chains),
#                                        # particles to be restraints
#                                        1.0,  # force constant
#                                        10.0,
#                                        # slack parameter affecting speed only
#                                        "Excluded-Volume"  # a string identifier
#                                        )
#
# """The springs between the beads of each protein chain were already defined by
# ProteinChainFactory, so now we create `rsf`, a scoring function
# (a.k.a. energy function) made of our restraints."""
#
# restraints = [chain.restraint for chain in chains] + [evr]
# rsf = IMP.core.RestraintsScoringFunction(restraints,
#                                          "Scoring function")  # Energy function
#
# """### III. Add dynamics
# Finally, we create a [`BrownianDynamics`]
# (https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1atom_1_1BrownianDynamics.html)
#  object, which enables us to simulate our interacting parts:
# """
#
# # BD
# bd = IMP.atom.BrownianDynamics(m)
# bd.set_scoring_function(rsf)  # tell BD about our scoring (energy) function
# bd_step_size_fs = 1000.0  # simulation time step in femotoseconds (10^-15 sec)
# bd.set_maximum_time_step(bd_step_size_fs)
# T = 300  # temperature in Kalvin
# bd.set_temperature(T)
#
# """#### Making a movie (trajectory) file
# We conclude by telling our simulation to output a trajectory (movie)
# file every 10000 frames. The movie is saved using the [Rich Molecular File format]
# (https://integrativemodeling.org/rmf/format.html), which can be viewed using e.g. the [Chimera software](https://www.cgl.ucsf.edu/chimera/download.html) or [ChimeraX](https://cxtoolshed.rbvi.ucsf.edu/apps/chimeraxrmf).
#
# **Note:** if you get a "UsageError: "Opening a file that is still being written
# is asking for trouble." error for any reson, just change the filename variable
# """
#
# rmf_filename = "my_trajectory.rmf"
# rmf = RMF.create_rmf_file(rmf_filename)
# rmf.set_description("Brownian dynamics trajectory with {}fs timestep.\n" \
#                     .format(bd_step_size_fs))
# IMP.rmf.add_hierarchy(rmf,
#                       h_root)  # Telling the movie that it should save all descendents of h_root
# IMP.rmf.add_restraints(rmf,
#                        restraints)  # the restraints are also saved and can be viewed in Chimera
# sos = IMP.rmf.SaveOptimizerState(m,
#                                  rmf)  # an optimizer state is invoked every n frames of simulation by the BD simulation
# sos.set_simulator(bd)
# frames_interval = 10000
# sos.set_period(frames_interval)
# bd.add_optimizer_state(sos)
# sos.update_always(
#     "initial conformation")  # save the initial conformation to the RMF file
#
# """## Simulating dynamics
# Without further ado, let's simulate our system using `bd.optimize(n)` where `n` is the number of time steps to simulate.
#
# We shall keep a few statistics every `n_inner_cycle` frames - the simulation time, the energy of the system (our scoring function), and the distance between the first and last bead of each chain.
#
# """
#
# T_ns = []  # time in nanoseconds
# E = []  # energy
# D = [[] for chain in chains]
# n_outer = 50000  # outer loop number of iterations
# n_inner = 250  # optimization per iteration
# for i in range(n_outer):
#     time_fs = bd.get_current_time()
#     time_ns = time_fs * 1e-6  # a nanosecond is a million femtoseconds
#     if i % (n_outer // 10) == 0:
#         print(f"Simulated for {time_ns:.1f} nanoseconds so far")
#     bd.optimize(n_inner)
#     T_ns.append(time_ns)  # keep time
#     E.append(bd.get_last_score())  # keep energy
#     for i, chain in enumerate(chains):  # keep distances for each chain
#         distance = IMP.core.get_distance(IMP.core.XYZ(chain.beads[0]),
#                                          IMP.core.XYZ(chain.beads[-1]))
#         D[i].append(distance)
# print(f"FINISHED. Simulated for {time_ns:.1f} nanoseconds in total.")
# rmf.flush()  # make sure RMF file is properly saved
#
# """Let's convert to numpy array - it will be easier to manipulate that way:"""
#
# T_ns = np.array(T_ns)  # this will be convenient later
# E = np.array(E)  # likewise
# D = np.array(D)
#
# """**Explanations:**
#
# To compute the energy of the system we use the `bd.get_last_score()` method.
#
# To compute the distance between the centers of two bead particles, we use
# the [`IMP.core.get_distance()`]
# (https://integrativemodeling.org/2.15.0/doc/ref/namespaceIMP_1_1core.html#af2049fdca8aae38e86549cfdeff1d1fc) function. This functions acts on [`IMP.core.XYZ`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1core_1_1XYZ.html) objects.
#
# *Advanced note:* bead particles can be converted (cast) to
# `IMP.core.XYZ`
# objects because they were decorated as such by our `ProteinChainFactory` class.
#
# ## Visualizing simulation results
# Let's see what happened. Since we kept statistics, we can now plot the energy and the distance between beads as a function of time. Optionally, you can frist inspect the resulting simulation as a movie in Chimera.
# """
#
# ### ADD CODE HERE #
# plt.xlabel(r'time [$ns$]')
# plt.ylabel(r'E [$kcal/mol/A^2$]')
# plt.show()
#
# # A3
#
# T_ns_threshold = 0.1  # CHANGE THIS LINE ONLY
# no_start = (T_ns > T_ns_threshold)
# plt.plot(T_ns[no_start], E[no_start], '-')
# plt.xlabel(r'time [$ns$]')
# plt.ylabel(r'E [$kcal/mol/A^2$]')
# plt.show()
#
# print(f"The mean energy is {np.mean(E[no_start]):.2f}")
# print(f"The std-deviation of the energy is {np.std(E[no_start]):.2f}")
#
# colors = cm.get_cmap("tab10")
# for i, Di in enumerate(D):
#     plt.plot(T_ns[no_start], Di[no_start],
#              color=colors(i % 10), label=f'D[{i}]')
# plt.xlabel(r'time [$ns$]')
# plt.ylabel(r'end-to-end distance [$A$]')
# plt.legend()
# plt.show()
#
# ax = sns.histplot(E[no_start])
# ax.set_xlabel(r'energy [$kcal^{-1}mol^{-1}A^2$]')
# plt.show()
#
# ax = sns.histplot(np.concatenate(D[0:]))
# ax.set_xlabel(r'N''-C'' distance [A]')
# plt.show()
#
# """####**Enrichment:**
# Let's inspect the code of classes `ProteinChain` and `ProteinChainFactory`
# way above by pressing "show code". This brief explanation is provided for
# context, you can move on with the exercise if you like.
#
# The method `create_bead()` creates a single bead, and decorates it with  [`IMP.core.XYZR`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1core_1_1XYZR.html), [`IMP.atom.Mass`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1atom_1_1Mass.html), [`IMP.atom.Diffusion`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1atom_1_1Diffusion.html), and [`IMP.display.Colored`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1display_1_1Colored.html).
#
#     def _create_bead(self,
#                      name : str):
#         p= IMP.Particle(self.model, name)
#         p_as_xyzr= IMP.core.XYZR.setup_particle(p) # A Decorator design pattern - adding functionality to an object at run time (~run-time inheritance)
#         p_as_xyzr.set_coordinates_are_optimized(True)
#         p_as_xyzr.set_radius(self.default_radius_A)
#         IMP.atom.Mass.setup_particle(p, FAKE_MASS)   # required by Hierarchy
#         IMP.atom.Diffusion.setup_particle(p)
#         IMP.atom.Hierarchy.setup_particle(p) # allow inclusion in IMP hierarchies
#         IMP.display.Colored.setup_particle(p,
#                                            IMP.display.get_display_color(0))
#         return p
#
#
# The method `_create_restraint` adds a spring to the simulation between each pair of consecutive beads. It first creates a function that increases quadratically from a certain rest distance - using the class [`IMP.core.HarmonicDistancePairScore`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1core_1_1HarmonicDistancePairScore.html). The [`IMP.container.PairsRestraint`](https://integrativemodeling.org/2.15.0/doc/ref/classIMP_1_1container_1_1PairsRestraint.html) object `pr` is a restraints that applies this function (`hdps`) over each pair of consecutive beads (`cpc`).
#
#     def _create_restraint(self,
#                           beads #spring constant
#                        ):
#       hdps = IMP.core.HarmonicDistancePairScore(self.rest_distance_A,
#                                                 self.k_kcal_per_mol_per_A2)
#       cpc = IMP.container.ConsecutivePairContainer(self.model, beads) # convention - use abbreviation for multiword class names
#       pr = IMP.container.PairsRestraint(hdps, cpc)
#       return pr, hdps
#
# # Probabilistic modeling (enrichment)
#
# There are no questions to answer in this section, so you can skip right on to **Part II** if you like. If you're still here, you will find some simple code for fitting the probability distribution of the energy and chain end-to-end distance distributions above, using the probabilistic modeling package [PyMC3](https://docs.pymc.io/).
#
# First install some required packages:
# """
#
# # !pip
# # install
# # theano - pymc
# # !pip
# # install
# # pymc3 == 3.11
# # .1
#
# """Define a probabilistic model made of a bunch of random variables ("RVs"):"""
#
# pm_model = pm.Model()
# with pm_model:
#     rv_E_mu = pm.Normal('E_mu', mu=20, sigma=20)
#     rv_E_sigma = pm.Normal('E_sigma', mu=20, sigma=20)
#     rv_E = pm.Normal('E', mu=rv_E_mu, sigma=rv_E_sigma,
#                      observed=E)
#     rv_D_mu = pm.Normal('D_mu', mu=100, sigma=50)
#     rv_D_sigma = pm.Normal('D_sigma', mu=50, sigma=25)
#     rv_D = pm.Normal('D', mu=rv_D_mu, sigma=rv_D_sigma,
#                      observed=np.concatenate(D))
# gv = pm.model_to_graphviz(pm_model)
#
# """Plot probabilistic relations:"""
#
# gv
#
# """Run Hamiltonian Markov-Chain Monte-Carlo to sample the mean (mu) and
# std-dev (sigma) for the various distributions:"""
#
# with pm_model:
#     trace = pm.sample(500, chains=2)
#
# # !pip install arviz
#
#
# with pm_model:
#     az.plot_trace(trace)
#
# az.summary(trace)
