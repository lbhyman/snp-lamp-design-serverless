import random as rnd
import json
from chalicelib.probe import Probe

# Randomly generate a starting population for the GA
def generate_initial_population(population_size, probe_params):
    print('Generating GA Population')
    probe_params['WT'] = probe_params['WT'].upper()
    probe_params['SNP'] = probe_params['SNP'].upper()
    probe_params['truncations'] = []
    probe_params['beta'] = [0,0,0,0]
    population = []
    for i in range(population_size):
        probe = Probe(probe_params)
        population.append(probe)
    return population

# Run a single generation of the GA
def generate_next_population(population):
    print('Generating GA Population')
    population.sort(reverse=True)
    if len(population) <= 2:
        return [population[0]]
    population_size = int(len(population)/2)
    num_parents = int(len(population)/4)
    # Ensure that all parents are unique
    parents = {}
    for probe in population:
        key = probe.get_key()
        if(len(parents) < num_parents):
            if(not key in parents):
                parents[key] = probe
    parents = list(parents.values())
    new_population = parents
    # Generate children
    for i in range(population_size - num_parents):
        parent_1 = rnd.choice(parents)
        parent_2 = rnd.choice(parents)
        # Improve diversity by introducing randomness when parents are identical
        if(parent_2.get_truncations() == parent_1.get_truncations()):
            random_probe_params = parent_1.get_dict()
            random_probe_params['truncations'] = []
            random_probe_params['beta'] = [0,0,0,0]
            parent_2 = Probe(random_probe_params)
            child = parent_1.cross(parent_2)
        else:
            child = parent_1.cross(parent_2)
        new_population.append(child)
    return new_population

# Calculate population fitness
def calculate_population_fitness(population):
    print('Calculating Fitness')
    for i in range(len(population)):
        curr_probe = population[i]
        if curr_probe.beta == [0,0,0,0]:
            curr_probe.calc_beta()
        population[i] = curr_probe
    population.sort(reverse=True)
    return population

# Generate next hill-climbing steps
def hill_climb(probe):
    print('Generating Hill-Climbing Population')
    return probe.next_iteration()