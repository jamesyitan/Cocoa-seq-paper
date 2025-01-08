import random
import numpy as np
import pandas as pd

def generate_unique_templates(num_unique_templates, num_bases):
    # List to store generated template molecules
    templates = []

    # Function to generate a random DNA sequence of the specified length
    def generate_random_sequence(length):
        bases = ['A', 'T', 'C', 'G']
        return ''.join(random.choice(bases) for _ in range(length))

    # Generate template molecules
    while len(templates) < num_unique_templates:
        template = generate_random_sequence(num_bases)
        templates.append(template)

    return templates

def generate_mock_community(species_list, molecules_per_species):
    # List to store generated mock community
    mock_community = []

    # Check if the lengths of species_list and molecules_per_species are the same
    if len(species_list) != len(molecules_per_species):
        raise ValueError("The lengths of species_list and molecules_per_species should be the same.")

    # Generate mock community
    for species, num_molecules in zip(species_list, molecules_per_species):
        # Add molecules for each species
        mock_community.extend([species] * num_molecules)

    # Shuffle the mock community to randomize the order
    random.shuffle(mock_community)

    return mock_community

def pcr(template_molecules, amplification_efficiency, num_cycles):
    # Initial pool of template molecules
    initial_pool = template_molecules.copy()

    # Simulation loop for each PCR cycle
    for cycle in range(num_cycles):
        # New pool for the current cycle
        current_pool = []

        # Iterate through each template molecule
        for template in initial_pool:
            # Determine amplification for the current template based on binomial distribution
            amplification = np.random.binomial(1, amplification_efficiency)

            # If amplification occurs, add the template to the current pool
            if amplification:
                current_pool.append(template)

        # Update the initial pool for the next cycle
        initial_pool = initial_pool + current_pool

    # Final pool of template molecules and amplicons
    final_pool = initial_pool

    return final_pool

def df_maker(pool_list, num_tracker):
    summary_df = pd.DataFrame(pool_list, columns=['Template'])
    occurrences_df = summary_df['Template'].value_counts().reset_index()
    occurrences_df.columns = ['Template', 'Occurrences']
    occurrences_df['Simulation'] = num_tracker

    return occurrences_df

def sim(num_sims, initial_templates, amplification_efficiency, num_cycles, seq_depth, save_tsv):
    tracker_df = pd.DataFrame()

    for sim in range(num_sims):
        final_pool = pcr(initial_templates, amplification_efficiency, num_cycles)
        subsamp_final_pool = random.sample(final_pool, seq_depth) #subsample to account for a standard seq depth
        final_df = df_maker(subsamp_final_pool,sim)
        tracker_df = pd.concat([tracker_df, final_df], ignore_index=True)
        tracker_df.to_csv(save_tsv, sep='\t', index=False)
    return final_pool

##RUNNING MORE THAN 20 CYCLES FOR PCR CAN LAST HOURS ON 200GB LARGE MEMORY NODE

##umi simulation
# initial_templates = generate_unique_templates(20, 6)
# sim_var_60eff = sim(5, initial_templates, 0.6, 30, 1000, 'sim_umi_seq/umi_20mol_60eff_30cycles_1000seq.tsv')
# sim_var_70eff = sim(10, initial_templates, 0.7, 30, 1000, 'sim_umi_seq/umi_20mol_70eff_30cycles_1000seq.tsv')
# sim_var_80eff = sim(5, initial_templates, 0.8, 30, 1000, 'sim_umi_seq/umi_20mol_80eff_30cycles_1000seq.tsv')
# sim_var_90eff = sim(5, initial_templates, 0.9, 30, 1000, 'sim_umi_seq/umi_20mol_90eff_30cycles_1000seq.tsv')
# sim_var_95eff = sim(5, initial_templates, 0.95, 30, 1000, 'sim_umi_seq/umi_20mol_95eff_30cycles_1000seq.tsv')

# ##community simulation
# #high biomass - uneven
# high_list = ['high-100', 'high-300', 'high-600']
# high_molecules_per_species = [100, 300, 600]
# high_initial = generate_mock_community(high_list, high_molecules_per_species)
# high_sim = sim(5, high_initial, 0.7, 20, 1000, 'sim_pcr_seq/high_uneven-100-300-600_70eff_20cycles_1000seq.tsv')

# #mid biomass - uneven
# mid_list = ['mid-10', 'mid-30', 'mid-60']
# mid_molecules_per_species = [10, 30, 60]
# mid_initial = generate_mock_community(mid_list, mid_molecules_per_species)
# mid_sim = sim(5, mid_initial, 0.7, 20, 1000, 'sim_pcr_seq/mid_uneven-10-30-60_70eff_20cycles_1000seq.tsv')

# #low biomass - uneven
# low_list = ['low-1', 'low-3', 'low-6']
# low_molecules_per_species = [1, 3, 6]
# low_initial = generate_mock_community(low_list, low_molecules_per_species)
# low_sim = sim(5, low_initial, 0.7, 20, 1000, 'sim_pcr_seq/low_uneven-1-3-6_70eff_20cycles_10000seq.tsv')

#extremely uneven but high biomass
uneven_high_list = ['uneven-1', 'uneven-99', 'uneven-900']
uneven_high_molecules_per_species = [1, 99, 900]
uneven_high_initial = generate_mock_community(uneven_high_list, uneven_high_molecules_per_species)
uneven_high_sim = sim(10, uneven_high_initial, 0.7, 30, 100000, 'sim_pcr_seq/uneven_high-1-99-900_70eff_30cycles_100000seq.tsv')
# #change code so you can do different subsampling schemes for the simulations

#extremely uneven but high biomass
# uneven_pair_list = ['uneven-1', 'uneven-100']
# uneven_pair_molecules_per_species = [1, 100]
# uneven_pair_initial = generate_mock_community(uneven_pair_list, uneven_pair_molecules_per_species)
# uneven_pair_sim_seq10 = sim(10, uneven_pair_initial, 0.7, 20, 10, 'sim_pcr_seq/uneven_pair-1-99_70eff_20cycles_10seq.tsv')
# uneven_pair_sim_seq100 = sim(10, uneven_pair_initial, 0.7, 20, 100, 'sim_pcr_seq/uneven_pair-1-99_70eff_20cycles_100seq.tsv')
# uneven_pair_sim_seq1000 = sim(10, uneven_pair_initial, 0.7, 20, 1000, 'sim_pcr_seq/uneven_pair-1-99_70eff_20cycles_1000seq.tsv')