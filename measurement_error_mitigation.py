# Default method:1 from [[file:~/Dropbox/Research/error_mitigation_nisq/scripts/org/measurement_error_mitigation.org::*Default method][Default method:1]]
default = 'ibu'
# Default method:1 ends here

# Other necessary modules:1 from [[file:~/Dropbox/Research/error_mitigation_nisq/scripts/org/measurement_error_mitigation.org::*Other necessary modules][Other necessary modules:1]]
import itertools
import numpy as np
def generate_all_bitstrings(num):
    """
    Generate all binary strings with num digits
    For n=2, this should return ['00', '01', '10', '11']
    """
    return [
        '0000'+''.join(x) for x in list(itertools.product(['0', '1'], repeat=num))
    ]
# Other necessary modules:1 ends here

# Response Matrix Method IBM:1 from [[file:~/Dropbox/Research/error_mitigation_nisq/scripts/org/measurement_error_mitigation.org::*Response Matrix Method IBM][Response Matrix Method IBM:1]]
def response_matrix_mitigation(M, y):
    M_array = np.array(M)
    M_inv = np.linalg.inv(M_array)
    y_array = np.array(y)
    y_mit = np.dot(M_inv, y_array)
    y_mit_list = y_mit.tolist()
    return y_mit_list
# Response Matrix Method IBM:1 ends here

# Vectorized error mitigation:1 from [[file:~/Dropbox/Research/error_mitigation_nisq/scripts/org/measurement_error_mitigation.org::*Vectorized error mitigation][Vectorized error mitigation:1]]
def unitary_mitigation(M, y):
    ytot = sum(
        y
    )  #for proper normalization that includes both counts and probabilites
    U = np.array([[np.sqrt(abs(Mij)) for Mij in Mj] for Mj in M])
    ys = np.array([np.sqrt(abs(yj)) for yj in y])
    Uinv = np.linalg.inv(U)
    xs = np.dot(Uinv, ys)
    norm = np.linalg.norm(xs)
    xs_hat = (xs / norm) * (np.sqrt(ytot))
    x = [xi**2 for xi in xs_hat.tolist()]
    return x
# Vectorized error mitigation:1 ends here

# Iterative Bayesian Bootstrapping:1 from [[file:~/Dropbox/Research/error_mitigation_nisq/scripts/org/measurement_error_mitigation.org::*Iterative Bayesian Bootstrapping][Iterative Bayesian Bootstrapping:1]]
def Rtj(Rmat, mvec, tvec, j):
    Rtj_vec = [Rmat[j][k] * tvec[k] for k in range(0, len(mvec))]
    return sum(Rtj_vec)


def Rtmi(Rmat, mvec, tvec, i):
    Rtmi_vec = [(Rmat[j][i] * tvec[i] * mvec[j]) / (Rtj(Rmat, mvec, tvec, j) + 10**(-15)) for j in range(0, len(mvec))]
    return sum(Rtmi_vec)


def tnvec(Rmat, mvec, tvec, n):
    for i in range(0, n):
        tnew = [Rtmi(Rmat, mvec, tvec, i) for i in range(0, len(mvec))]
        tvec = tnew
    return tvec


def iterative_bayesian_mitigation(M, y):
    """
    M is the response matrix
    y are the observed counts vector (not probability vector)
    """
    y_unit = [x / sum(y) for x in y] #probability vector
    y_mit = tnvec(M, y_unit, y_unit, 10) #mitigated probability vector
    return [x*sum(y) for x in y_mit] #mitigated counts vector
# Iterative Bayesian Bootstrapping:1 ends here

# Get MEM matrix from data:1 from [[file:~/Dropbox/Research/error_mitigation_nisq/scripts/org/measurement_error_mitigation.org::*Get MEM matrix from data][Get MEM matrix from data:1]]
def mem_matrix_from_data(mem_data,num_qubits=1):
    """
    mem_data is a list of dictionaries. For n qubits, this should have 2^n different datasets.
    The list should be organized in increasing bitstring order - like 00, 01, 10, 11
    Output: a matrix which is simply a list of lists like [[1,2], [3,4]]
    
    'matrix' defined below gives C_noisy = transpose(matrix) C_ideal 
    the mem_matrix = inverse(transpose(matrix)) gives C_ideal = mem_matrix C_noisy
    """
    base_data = mem_data[0]  #first dataset in the list
    base_string = list(base_data.keys())[0]  #first bitstring in basedata
    # num_qubits = len(base_string)
    assert (len(mem_data) == 2**num_qubits)
    all_strings = generate_all_bitstrings(num_qubits)
    matrix = []
    total = sum(list(
        base_data.values()))  #use the total to go from counts to probabilities
    for data in mem_data:  #iterate over the input strings
        row = []
        for string in all_strings:
            value = data[string]
            row.append(value)
        matrix.append(row)
    matrix = np.array(matrix) / total
    mem_matrix = np.transpose(matrix)
    # mem_matrix = np.linalg.inv(matrix)
    return mem_matrix
# Get MEM matrix from data:1 ends here

# Get MEM matrix from data:2 from [[file:~/Dropbox/Research/error_mitigation_nisq/scripts/org/measurement_error_mitigation.org::*Get MEM matrix from data][Get MEM matrix from data:2]]
def mem_matrix_from_data_test():
    mem_data = [{
        '10': 96,
        '11': 1,
        '01': 95,
        '00': 9808
    }, {
        '10': 2,
        '11': 103,
        '01': 9788,
        '00': 107
    }, {
        '10': 9814,
        '11': 90,
        '01': 1,
        '00': 95
    }, {
        '10': 87,
        '11': 9805,
        '01': 107,
        '00': 1
    }]
    return mem_matrix_from_data(mem_data)


# mem_matrix_from_data_test()
# Get MEM matrix from data:2 ends here

# Choosing the appropriate method:1 from [[file:~/Dropbox/Research/error_mitigation_nisq/scripts/org/measurement_error_mitigation.org::*Choosing the appropriate method][Choosing the appropriate method:1]]
def mem_on_data(data, mem_matrix, method = default,num_qubits=1):
    """
    data : dictionary
    mem_matrix : a list of lists (not array)
    return data dictionary with adjusted counts
    """
    method_dic = {'vec':unitary_mitigation, 'ibu':iterative_bayesian_mitigation, 'ibm':response_matrix_mitigation}

    # num_qubits = len(list(data.keys())[0])
    bitstrings = generate_all_bitstrings(num_qubits)
    data_vector = []  # of the form [a, b, c, d]
    
    for string in bitstrings:
        data_vector.append(data[string])

    # print(method)
    if method == 'vec':
        corrected_data_list = unitary_mitigation(mem_matrix, data_vector)

    if method == 'ibm':
        corrected_data_list = response_matrix_mitigation(mem_matrix, data_vector)

    if method == 'ibu':
        corrected_data_list = iterative_bayesian_mitigation(mem_matrix, data_vector)
    
    new_dict = {}
    for index in range(0, len(bitstrings)):
        string = bitstrings[index]
        corrected_count = corrected_data_list[index]
        new_dict[string] = corrected_count
    return new_dict
# Choosing the appropriate method:1 ends here
