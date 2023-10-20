# TODO: "THE GAP"
# new model, set Kc = 0 and set other K to 1, change Kf accordingly
# change end time, i.e. t = 10, 100, ...
# Ki is constant, others are linear

import numpy as np
import random
import math
from collections import defaultdict

def calculate_intensity(reaction_count, reactants_list, species_counts, rates, species_list):
    intensities = []
    for i in range(reaction_count):
        rate = rates[i]
        for reactant in reactants_list[i]:
            species, coefficient = reactant
            
            if coefficient == 1:
                index = species_list.index(species)
                rate *= species_counts[index]
            elif coefficient == 2:
                index = species_list.index(species)
                rate *= 0.5 * species_counts[index] * (species_counts[index] - 1)
            # elif coefficient == 0:
            #     index = species_list.index(species)
            #     rate *= rates[index]
        intensities.append(rate)
    return intensities

def reaction_to_changes(reactants, products):
    changes = defaultdict(int)  # Initialize a dictionary to store changes for each species
    
    # Process reactants
    for reactant in reactants:
        species, coefficient = reactant
        changes[species] -= coefficient
    
    # Process products
    for product in products:
        species, coefficient = product
        changes[species] += coefficient
    
    # Convert the dictionary to a list
    changes_list = [(species, change) for species, change in changes.items()]
    result = [x[1] for x in changes_list]
    print(f'The jumps of reaction: {result}')
    return result

def simGillespie(species_counts, jumps, rates, end_time, reaction_count, reactants_list, species_list):
    # Define the initial species counts, reaction rates, and intensity functions
    # species_counts = [10, 10, 50, 10]  # Initial counts for four species
    # rates = [200, 10, 25, 1, 0.01, 1]  # Reaction rates for six reactions

    # Define the changes in species counts for each reaction
    # jumps = [
    #     [0, 1, 0, 0],      # Transcription
    #     [0, 0, 1, 0],      # Translation
    #     [0, -1, 0, 0],     # Degradation of mRNA
    #     [0, 0, -1, 0],     # Degradation of protein
    #     [0, 0, -2, 1],     # Protein dimerization
    #     [0, 0, 0, -1],     # Degradation of dimer
    # ]

    # Define the time and simulation end time
    time = 0
    break_early = 0
    # end_time = 8

    # Initialize lists to store time points and species counts
    time_points = [time]
    species_counts_history = [species_counts]

    # Main simulation loop
    while time < end_time:
        intensity_functions = calculate_intensity(reaction_count, reactants_list, species_counts, rates, species_list)
        # print(f'Intensities function: {intensity_functions}')
        # Calculate the intensity functions based on species_counts
        # intensity_functions = [rates[0], 
        #                     rates[1] * species_counts[1], 
        #                     rates[2] * species_counts[1], 
        #                     rates[3] * species_counts[2], 
        #                     rates[4] * species_counts[2] * (species_counts[2] - 1), 
        #                     rates[5] * species_counts[3]]
        
        # Calculate the total intensity
        total_intensity = sum(intensity_functions)
        
        # Generate two independent uniform(0,1) random numbers
        r1 = random.uniform(0, 1)
        r2 = random.uniform(0, 1)
        
        # Calculate the time increment
        delta_t = -np.log(r1) / total_intensity
        
        # TODO: what to do if total_intensity is 0 in Gillespie's algorithm?
        if delta_t <= 0 or total_intensity == 0:
            break_early = 1
            break
        # print(f'Total_intensity is : {total_intensity}, Delta_t for reaction is: {delta_t}')
        
        # Find index
        i = 0
        p_sum = 0.0
        while p_sum < r2 * total_intensity:
            p_sum += intensity_functions[i]
            i += 1
        reaction = i - 1
        
        # Update time
        time += delta_t
        
        # Update species counts based on the chosen reaction
        species_counts = [count + change for count, change in zip(species_counts, jumps[reaction])]
        
        # Append the current time and species counts to the history
        time_points.append(time)
        species_counts_history.append(list(species_counts))
    
    # return species_counts
    if break_early == 1:
        return species_counts
    else:
        return species_counts_history[len(species_counts_history) - 2]



def full_model(init_comp_count, num_species, init_species_counts, network_jumps, network_rates, end_time, reaction_count, reactants_list, species_list):
    
    comp_counts = init_comp_count
    comp_rates = [1, 1, 1, 1]  # Reaction rates for four compartment reactions
    curr_state = init_species_counts

    # Define the changes in species counts for each reaction
    # Reaction 0: 0 -> C
    # Reaction 1: C -> 0
    # Reaction 2: 2C -> C
    # Reaction 3: C -> 2C
    comp_jumps = [1, -1, -1, 1]

    # Define the time and simulation end time
    time = 0
    end_time = 100

    # Initialize lists to store time points and species counts
    time_points = [time]
    # species_counts_history = [comp_counts]

    # Main simulation loop
    while time < end_time:
        # Calculate the intensity functions based on species_counts
        intensity_functions = [comp_rates[0], 
                            comp_rates[1] * comp_counts, 
                            comp_rates[2] * comp_counts * (comp_counts - 1) / 2,
                            comp_rates[3] * comp_counts]
        
        # Calculate the total intensity
        total_intensity = sum(intensity_functions)
        
        # Generate two independent uniform(0,1) random numbers
        r1 = random.uniform(0, 1)
        r2 = random.uniform(0, 1)
        
        # Calculate the time increment
        delta_t = -np.log(r1) / total_intensity
        
        # Update time
        time += delta_t
        
        # Decide which reaction occurred
        i = 0
        p_sum = 0.0
        while p_sum < r2 * total_intensity:
            p_sum += intensity_functions[i]
            i += 1
        reaction = i - 1

        # Update species counts based on the chosen reaction
        comp_counts = comp_counts + comp_jumps[reaction]

        # Append the current time and species counts to the history
        time_points.append(time)
        # species_counts_history.append(comp_counts)

        # 1. if 0 -> C is chosen
        if reaction == 0:
            print("0 -> C is chosen")
            for i in range(comp_counts - 1):
                newState = simGillespie(curr_state[i], network_jumps, network_rates, delta_t, reaction_count, reactants_list, species_list)
                curr_state[i] = newState

            # the initial state of the new compartment
            # TODO: should be determined by miu
            new_comp = [random.randint(0, 10) for _ in range(num_species)]
            curr_state.append(new_comp)
        
        # 2. if C -> 0 is chosen
        elif reaction == 1:
            print("C -> 0 is chosen")
            del_comp = random.randint(0, comp_counts)
            curr_state.pop(del_comp)

            for i in range(comp_counts):
                newState = simGillespie(curr_state[i], network_jumps, network_rates, delta_t, reaction_count, reactants_list, species_list)
                curr_state[i] = newState
        
        # 3. if 2C -> C is chosen
        elif reaction == 2:
            print("2C -> C is chosen")
            for i in range(comp_counts + 1):
                newState = simGillespie(curr_state[i], network_jumps, network_rates, delta_t, reaction_count, reactants_list, species_list)
                curr_state[i] = newState

            # choose two different compartments uniformly at random
            two_random_comp = []

            while len(two_random_comp) < 2:
                random_int = random.randint(0, comp_counts)
                if random_int not in two_random_comp:
                    two_random_comp.append(random_int)

            curr_state[min(two_random_comp)] = [x + y for x, y in zip(curr_state[two_random_comp[0]], curr_state[two_random_comp[1]])]
            curr_state.pop(max(two_random_comp))
        
        # 4. if C -> 2C is chosen
        else:
            print("C -> 2C is chosen")
            for i in range(comp_counts - 1):
                newState = simGillespie(curr_state[i], network_jumps, network_rates, delta_t, reaction_count, reactants_list, species_list)
                curr_state[i] = newState

            # choose a compartment to split uniformly at random
            split_comp = random.randint(0, comp_counts - 2)

            new_comp_1 = list(curr_state[split_comp])
            new_comp_2 = list(curr_state[split_comp])

            for j in range(len(curr_state[split_comp])):
                new_comp_1[j] = random.randint(0, curr_state[split_comp][j])
                new_comp_2[j] = new_comp_2[j] - new_comp_1[j]
            
            curr_state.pop(split_comp)
            curr_state.append(new_comp_1)
            curr_state.append(new_comp_2)

        print(f'Current state is {curr_state}')
    
    return curr_state


def main():
    # Ask for user input and parse it
    init_comp_count = int(input("Enter the initial compartment count: "))
    species_input = input("Enter all species in the reaction network (e.g., 'G M P'): ")
    species_list = species_input.split()
    num_species = len(species_list)
    # num_species = int(input("Enter the number of species in the reaction network: "))

    # Parse a list of initial species counts
    if init_comp_count == 0:
        init_species_counts = []
    else:
        init_species_counts = []
        for i in range(init_comp_count):
            input_species_counts = input(f"Enter initial species counts of compartment {i + 1} separated by spaces: ")
            init_species_counts.append([int(x) for x in input_species_counts.split()])
        
    reaction_count = int(input("Enter the number of reactions in reaction network: "))

    # Parse reaction jumps
    network_jumps = []
    reactants_list = []
    for i in range(reaction_count):
        # Ask for user input for reactants
        reactants = []
        reactants_input = input("Enter reactants separated by '+' (Ex. 2 A + 1 B): ")
        if reactants_input:
            reactant_tokens = reactants_input.split('+')
            for token in reactant_tokens:
                coefficient, element = token.strip().split()
                reactants.append((element, int(coefficient)))
        reactants_list.append(reactants)

        # Ask for user input for products
        products = []
        products_input = input("Enter products separated by '+' (Ex. 2 A + 1 B): ")
        if products_input:
            product_tokens = products_input.split('+')
            for token in product_tokens:
                coefficient, element = token.strip().split()
                products.append((element, int(coefficient)))
        
        network_jumps.append(reaction_to_changes(reactants, products))

    # # Parse a list of network jumps
    # network_jumps = []
    # for i in range(reaction_count):
    #     input_jumps = input(f"Enter jump of reaction {i + 1} separated by spaces: ")
    #     network_jumps.append([int(x) for x in input_jumps.split()])

    # Parse a list of network rates
    network_rates = input("Enter network rates separated by spaces: ")
    network_rates = [float(x) for x in network_rates.split()]

    end_time = float(input("Enter the end time: "))

    output = full_model(init_comp_count, num_species, init_species_counts, network_jumps, network_rates, end_time, reaction_count, reactants_list, species_list)
    print(f"Final state is {output}")

if __name__ == '__main__':
    main()