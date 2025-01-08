"""Create sets for the optimization. The sets are used to define the variables
and constraints.
"""

import numpy as np


def create(matrices):
    """Create sets for the optimization. The sets are used to define the
    variables and constraints. Depending on the matrices, the sets are
    calculated with the defined directions (i -> j as a reference).

    Args:
        matrices (dict): Dictionary with the matrices of the district heating
        network with keys:
            * a_i: Incidence matrix with rows: nodes, columns: edges.
            * a_p: Adjacency matrix for the producers: rows: nodes, columns: supply ids.
            * a_c: Adjacency matrix for the consumers: rows: nodes, columns: consumer ids.

    Returns:
        dict: sets for the optimization.
            * connection_c_ij: Consumers connected to the network in direction i -> j
            * lambda_c_ij: Binary vector with 1 where consumers are connected
                to the network in direction i -> j
            * connection_c_ji: Consumers connected to the network in direction
                j -> i
            * lambda_c_ji: Binary vector with 1 where consumers are connected
                to the network in direction j -> i
            * a_i_out: Dictionary with the outgoing pipes for each node (
                direction i -> j is 1)
            * a_i_in: Dictionary with the ingoing pipes for each node
                (direction i -> j ingoing is -1)
            * a_p_in: Dictionary with the ingoing pipes for each producer
                (direction i -> j ingoing is -1)
            * a_c_out: Dictionary with the outgoing pipes for each consumer
                (direction i -> j is 1)
            * a_c_out_edge: Dictionary with the outgoing pipes for each consumer
                with the edge pipes
    """
    s = {}
    # Shapes of matrices
    s['a_i_shape'] = np.shape(matrices['a_i'])  # (rows 0, columns 1)
    s['a_p_shape'] = np.shape(matrices['a_p'])
    s['a_c_shape'] = np.shape(matrices['a_c'])

    # where consumers are connected to the network
    consumers = np.where(matrices['a_c'].sum(axis=1) == 1)[0]

    s['connection_c_ij'] = np.where(
        matrices['a_i'][consumers, :].sum(axis=0) == -1
    )[0]
    s['lambda_c_ij'] = np.zeros(s['a_i_shape'][1])
    s['lambda_c_ij'][s['connection_c_ij']] = 1

    s['connection_c_ji'] = np.where(
        matrices['a_i'][consumers, :].sum(axis=0) == 1
    )[0]
    s['lambda_c_ji'] = np.zeros(s['a_i_shape'][1])
    s['lambda_c_ji'][s['connection_c_ji']] = 1

    # create a dict with {row: [outgoing pipes], row: [ingoing pipes]}
    s['a_i_out'] = {}
    s['a_i_in'] = {}
    for i in range(s['a_i_shape'][0]):
        s['a_i_out'][i] = np.where(matrices['a_i'][i, :] == 1)[0]
        s['a_i_in'][i] = np.where(matrices['a_i'][i, :] == -1)[0]

    s['a_p_in'] = {}
    for i in range(s['a_p_shape'][0]):
        s['a_p_in'][i] = np.where(matrices['a_p'][i, :] == -1)[0]

    s['a_c_out'] = {}
    s['a_c_out_edge'] = {}
    for i in range(s['a_c_shape'][0]):
        s['a_c_out'][i] = np.where(matrices['a_c'][i, :] == 1)[0]
        if matrices['a_c'][i, :].any() == 1:
            s['a_c_out_edge'][i] = np.where(
                (matrices['a_i'][i, :] == 1)
                | (matrices['a_i'][i, :] == -1)
                )[0]
        else:
            s['a_c_out_edge'][i] = []

    return s
