import qutip
import random 
import numpy as np
from copy import deepcopy 


def xeb(rho_1, rho_2):
    return rho_1.shape[0]*(np.diagonal(rho_1).dot(np.diagonal(rho_2))) - 1

def xeb_log(rho_1, rho_2):
    return -1*rho_1.shape[0]*(np.diagonal(rho_1).dot(np.log(np.real(np.diagonal(rho_2)))))


def get_all_to_all_random_haar_circuit(N:int, tau:int, nearest_neigbhor=False):
    all_pairs = [(i, j) for i in range(N) for j in range(i)]
    gate_s = []
    for t in range(tau):
        if nearest_neigbhor is False:
            for (i, j) in random.sample(all_pairs, N//2):
                gate_s.append((i, j, random.randint(0, 100000), t))
        else:
            for i in range(t % 2, N, 2):
                qubit_a, qubit_b =i,(i+1)%N
                gate_s.append((qubit_a, qubit_b, random.randint(0, 10000), t))
    return gate_s



def two_ancillas_opposite_end(N, circuit, sigma=0.0, quenched=False, noise_samples=1):
    # in this case, apply the error and the inverse channel on the same copy 
    # sample each circuit for many samples 
    # return only the density matrix corresponding to the two qubits 

    
    alpha = sigma  # ratio (1-q1)/(1-qbar)
    qmean = 0.1  # average of q1/q2
    p = 0.9
    
    

    q1 = qmean - sigma*np.sqrt((1-p)/p) # such that p*q1 + (1-p)*q2 = qmean
    q2 = (qmean-p*q1)/(1-p)           

    qbar = 1-(np.abs(1-q1)**p)*(np.abs(1-q2)**(1-p))

 
    psi:qutip.Qobj = qutip.qip.qubit_states(N, [0]*(N))

    # entangle qubits at opposite ends 
    # apply a random Haar unitary circuit to get a completely scrambled state in the beginning 
    for t in range(N):
        for k in range((N)//2):
            i,j = np.random.choice(np.arange(N), replace=False, size=2)
            # apply random gate
            psi = qutip.operations.expand_operator(qutip.rand_unitary(4, dims=[[2, 2], [2, 2]]), N, [i, j])*psi 

    


    # helper objects 
    op_dict = {}
    for i in range(N):
        op_dict[i]=[]
        for pauli_op in [qutip.sigmax(), qutip.sigmay(), qutip.sigmaz()]:
            op_dict[i].append(qutip.operations.expand_operator(pauli_op, N, [i]))

    def depolarization_channel(rho, k, p):
        rho_tr = np.trace(rho)
        rho = (1-p + rho_tr*p/4)*rho + rho_tr*(p/4)*sum([op_dict[k][pauli_idx] * rho * op_dict[k][pauli_idx] for pauli_idx in range(3)])
        return rho 

    def inverse_depolarization_channel(rho, k, q1bar, q2bar):
        rho_tr = np.trace(rho)
        rho = (1/(1-q1bar))*((1-q2bar*rho_tr/4)*rho - rho_tr*(q2bar/4)*sum([op_dict[k][pauli_idx] * rho * op_dict[k][pauli_idx] for pauli_idx in range(3)]))
        return rho 

    if quenched:
        noise_dict = [[q1, q2][np.random.binomial(1, p)] for x in range(N)]


    # Now start the noisy dynamics. 
    rho_list = []

    for _ in range(noise_samples):
        rho = qutip.ket2dm(deepcopy(psi))
        for gate in circuit:
            i,j , seed, t = gate
            # there are two entangled qubits at either end, shift the indices by 1 by 1.

            # random unitary gate 
            unitary = qutip.operations.expand_operator(qutip.rand_unitary(4, seed=seed, dims=[[2, 2], [2, 2]]), N, [i, j])
            
            #psi = unitary * psi
            rho = unitary * rho * unitary.dag()

            for k in [i, j]:
                # apply depolarizing channel to each of the two qubits
                if not quenched:
                    qit = [q1, q2][np.random.binomial(1, 1-p)] # choose between q1/q2
                    # q1 with probability p and q2 with probabilty 1-p

                else:
                    qit = noise_dict[k]

                # noise channel
                rho = depolarization_channel(rho, k, qit)
                # antinoise channel 
                rho = rho/np.trace(rho)
                rho = inverse_depolarization_channel(rho, k, qbar, qmean)
                rho = rho/np.trace(rho)

        rho_list.append(rho.ptrace([0, N-1]))

    return rho_list 

