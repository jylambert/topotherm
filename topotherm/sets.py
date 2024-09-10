"""Create sets for the optimization. The sets are used to define the variables
and constraints.
"""
import numpy as np

def create(matrices):
    """Create sets for the optimization. The sets are used to define the
    variables and constraints.

    Args:
        matrices (dict): Dictionary with the matrices of the district heating
        network with keys a_i, a_p, a_c.

    Returns:
        s: dictionary with the sets for the optimization
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
    for i in range(s['a_c_shape'][0]):
        s['a_c_out'][i] = np.where(matrices['a_c'][i, :] == 1)[0]
    return s
