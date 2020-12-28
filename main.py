# usr/bin/env python3
"""runs different parts of the program."""
import src.performance as perf
import src.run_simulation as rs


# perf.particle_build_up(OM=8)

for dark in [False, True]:
    rs.run_through_galaxy_sizes(halo=dark)

# rs.run_through_impactors()

