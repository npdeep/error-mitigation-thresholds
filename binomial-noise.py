import utils 
from tqdm import tqdm 
import pickle 
from multiprocessing import Pool 
import numpy as np 

noise_samples=20
n_instance = 100
quenched = False 
nearest_neigbhor=quenched

def get_data(disorder_and_L):
    L, disorder = disorder_and_L
    xeb_data_for_disorder =[]
    for circuit_idx in tqdm(range(n_instance)):
        circuit = utils.get_all_to_all_random_haar_circuit(L, 4*L, nearest_neigbhor=nearest_neigbhor)    
        rho_noisy_s = utils.two_ancillas_opposite_end(L, circuit,  sigma=disorder, quenched=quenched, noise_samples=noise_samples)
        xeb_data_for_disorder.append(rho_noisy_s)
    return xeb_data_for_disorder 

def main():
    xeb_data = []
    disorder_s = np.linspace(0.05, 1, 23)

    for L in [4, 6, 8]:
        with Pool(23) as p:
            xeb_data_all = p.map(get_data, [(L, disorder) for disorder in disorder_s])

            # for circuit_idx in tqdm(range(n_instance)):
            #     circuit = utils.get_all_to_all_random_haar_circuit(L, 4*L, nearest_neigbhor=nearest_neigbhor)    
            #     rho_noisy_s = utils.two_ancillas_opposite_end(L, circuit, 
            #     disorder=disorder, quenched=quenched, noise_samples=noise_samples)
            #     xeb_data_for_disorder.append(rho_noisy_s)
            # xeb_data_all.append(xeb_data_for_disorder)

        with open(f"data/binomial-disorder-all-to-all-L={L}.pickle", "wb") as file:
            pickle.dump(xeb_data_all, file)

if __name__ == "__main__":
    main()