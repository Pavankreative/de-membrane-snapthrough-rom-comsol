"""Regenerate manuscript figures using the revision 17 archived inputs."""
from make_rev17_figures import trajectories, basins, manifolds, sections

if __name__ == '__main__':
    trajectories()
    basins()
    manifolds()
    sections()
