# Significance clustering


## Author:

Martin Rosvall

For contact information, see http://www.mapequation.org/about


## Getting started:

In a terminal with the GNU Compiler Collection installed,
just run 'make' in the current directory to compile the
code with the included Makefile.

Call: ./dangling-lumping [-s <seed>] input_state_network.net output_state_network.net
seed: Any positive integer.
input_state_network.net: The state network with state nodes with dangling state nodes (no out-links)
output_state_network.net: The lumped state network where all state nodes have been randomly merged
                          with non-dangling state nodes of the same physical node. If no such nodes
                          exist, the dangling nodes are lumped into a single dangling state node
                          per physics node.