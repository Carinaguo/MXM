# TODO: "THE GAP"
# new model, set Kc = 0 and set other K to 1, change Kf accordingly (2.5)
# change end time, i.e. t = 10, 100, ...
# Ki is constant, others are linear
# initial: 3 compartments with 1 S inside
# miu: random.randint(0,10)
# end time: 5

import numpy as np
import random
from collections import defaultdict
import matplotlib.pyplot as plt

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


def full_model(init_comp_count, init_species_counts, end_time):
    
    history = [0]
    comp_counts = init_comp_count
    reaction_rates = [1, 1, 0, 2.5, 1, 1] # ki, ke, kc, kf, kb, kd
    curr_state = init_species_counts

    S_counts = 0 # Total number of species S in all compartments
    for i in range(comp_counts):
        S_counts += curr_state[i][0]

    # Define the changes in species counts for each reaction
    # Reaction 0: 0 -> C
    # Reaction 1: C -> 0
    # Reaction 2: 2C -> C
    # Reaction 3: C -> 2C
    # Reaction 4: 0 -> S
    # Reaction 5: S -> 0
    comp_jumps = [1, -1, -1, 1]

    # Define the time and simulation end time
    time = 0
    # end_time = 100

    # Initialize lists to store time points and species counts
    time_points = [time]
    # species_counts_history = [comp_counts]

    # Main simulation loop
    while time < end_time:
        # Calculate the intensity functions based on species_counts
        intensity_functions = [reaction_rates[0], 
                            reaction_rates[1] * comp_counts, 
                            0,
                            reaction_rates[3] * S_counts,
                            reaction_rates[4] * comp_counts,
                            reaction_rates[5] * S_counts]
        
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
        if reaction != 4 and reaction != 5:
            comp_counts = comp_counts + comp_jumps[reaction]

        # Append the current time and species counts to the history
        time_points.append(time)
        # species_counts_history.append(comp_counts)

        # 1. if 0 -> C is chosen
        if reaction == 0:
            print("0 -> C is chosen")

            # the initial state of the new compartment
            # TODO: should be determined by miu
            new_comp = [random.randint(0, 10)]
            curr_state.append(new_comp)
            S_counts += new_comp[0]
        
        # 2. if C -> 0 is chosen
        elif reaction == 1:
            print("C -> 0 is chosen")
            del_comp = random.randint(0, comp_counts)
            S_counts -= curr_state[del_comp][0]
            curr_state.pop(del_comp)
        
        # 3. if 2C -> C is chosen
        elif reaction == 2:
            print("2C -> C is chosen")

            # choose two different compartments uniformly at random
            two_random_comp = []

            # if the number of comparten

            while len(two_random_comp) < 2:
                random_int = random.randint(0, comp_counts)
                if random_int not in two_random_comp:
                    two_random_comp.append(random_int)

            curr_state[min(two_random_comp)] = [x + y for x, y in zip(curr_state[two_random_comp[0]], curr_state[two_random_comp[1]])]
            curr_state.pop(max(two_random_comp))
            # S_counts stays the same
        
        # 4. if C -> 2C is chosen
        elif reaction == 3:
            print("C -> 2C is chosen")
            # for i in range(comp_counts - 1):
            #     newState = simGillespie(curr_state[i], network_jumps, network_rates, delta_t, reaction_count, reactants_list, species_list)
            #     curr_state[i] = newState

            # choose a compartment to split uniformly at random
            # split_comp = random.randint(0, comp_counts - 2)
            split_comp = np.random.choice(np.arange(comp_counts - 1), p = [curr_state[i][0] / S_counts for i in range(comp_counts - 1)]) # TODO

            new_comp_1 = list(curr_state[split_comp])
            new_comp_2 = list(curr_state[split_comp])

            for j in range(len(curr_state[split_comp])):
                new_comp_1[j] = random.randint(0, curr_state[split_comp][j])
                new_comp_2[j] = new_comp_2[j] - new_comp_1[j]
            
            curr_state.pop(split_comp)
            curr_state.append(new_comp_1)
            curr_state.append(new_comp_2)
            # S_counts stays the same

        # 5. if 0 -> S is chosen
        elif reaction == 4:
            print("0 -> S is chosen")
            comp = random.randint(0, comp_counts-1)
            curr_state[comp][0] += 1
            S_counts += 1

        # 6. if S -> 0 is chosen
        elif reaction == 5:
            print("S -> 0 is chosen")
            comp = np.random.choice(np.arange(comp_counts), p = [curr_state[i][0] / S_counts for i in range(comp_counts)]) # TODO
            curr_state[comp][0] -= 1
            S_counts -= 1


        print(f'Current state is {curr_state}')
        history.append(S_counts)
    
    print(f'Time points: {time_points}')
    # Create a basic line plot
    plt.plot(time_points, history)

    # Add labels and a title
    plt.xlabel('time')
    plt.ylabel('S_counts')

    # Display the plot
    plt.show()
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
        
    # reaction_count = int(input("Enter the number of reactions in reaction network: "))

    # # Parse reaction jumps
    # network_jumps = []
    # reactants_list = []
    # for i in range(reaction_count):
    #     # Ask for user input for reactants
    #     reactants = []
    #     reactants_input = input("Enter reactants separated by '+' (Ex. 2 A + 1 B): ")
    #     if reactants_input:
    #         reactant_tokens = reactants_input.split('+')
    #         for token in reactant_tokens:
    #             coefficient, element = token.strip().split()
    #             reactants.append((element, int(coefficient)))
    #     reactants_list.append(reactants)

    #     # Ask for user input for products
    #     products = []
    #     products_input = input("Enter products separated by '+' (Ex. 2 A + 1 B): ")
    #     if products_input:
    #         product_tokens = products_input.split('+')
    #         for token in product_tokens:
    #             coefficient, element = token.strip().split()
    #             products.append((element, int(coefficient)))
        
    #     network_jumps.append(reaction_to_changes(reactants, products))

    # # Parse a list of network jumps
    # network_jumps = []
    # for i in range(reaction_count):
    #     input_jumps = input(f"Enter jump of reaction {i + 1} separated by spaces: ")
    #     network_jumps.append([int(x) for x in input_jumps.split()])

    # Parse a list of network rates
    # network_rates = input("Enter network rates separated by spaces: ")
    # network_rates = [float(x) for x in network_rates.split()]

    end_time = float(input("Enter the end time: "))

    output = full_model(init_comp_count, init_species_counts, end_time)
    print(f"Final state is {output}")

if __name__ == '__main__':
    main()