import random as rnd
import json
from probe import Probe

# Randomly generate a starting population for the GA
def generate_initial_population(self, population_size, probe_params):
    probe_params['truncations'] = []
    probe_params['beta'] = [0,0,0,0]
    population = []
    for i in range(population_size):
        parent = Probe(probe_params)
        population.append(parent.__dict__)
    return population

# Run a single generation of the GA
def generate_next_population(self, population):
    population = [Probe(probe_params) for probe_params in population]
    population.sort(reverse=True)
    if len(population) <= 2:
        return population
    population_size = len(population)/2
    num_parents = len(population)/4
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
            random_probe_params = parent_1.__dict__
            random_probe_params['truncations'] = []
            random_probe_params['beta'] = [0,0,0,0]
            parent_2 = Probe(random_probe_params)
        else:
            child = parent_1.cross(parent_2)
        new_population.append(child)
    return [p.__dict__ for p in new_population]