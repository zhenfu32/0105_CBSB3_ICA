import numpy as np

def selex_sequence_alignment(data):
    # Inputs:
    # data - a numpy array of sequences from all rounds of SELEX experiments
    # Returns:
    # aligned_data - a numpy array of aligned sequences from all rounds of SELEX experiments

    # Determine length of the binding site and the number of sequences
    l = len(data[0]) # length of the binding site
    N = data.shape[0] # number of sequences

    # Initialize an array for aligned sequences
    aligned_data = np.zeros((N,l))

    # Implement global sequence alignment
    for i in range(N):
        # Extract sequence i and reverse complement if necessary
        seq_i = data[i,:]
        if np.sum(seq_i == 0) <= l/2:
            seq_i = revcomp(seq_i) # reverse complement sequence if necessary

        # Perform global sequence alignment with respect to first round consensus sequence
        aligned_i, cost_i = global_align(seq_i,data[0,:])

        # Store aligned sequence
        aligned_data[i,:] = aligned_i[0]

    return aligned_data.astype(int)

def revcomp(seq):
    # Inputs:
    # seq - a numpy array of nucleotide sequence
    # Returns:
    # a numpy array of the reverse complement of the input sequence

    # Define complement dictionary
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # Reverse and complement sequence
    rev_seq = seq[::-1]
    revcomp_seq = np.array([comp_dict[nt] for nt in rev_seq])

    return revcomp_seq

