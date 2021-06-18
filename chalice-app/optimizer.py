from probe import Probe
import ga_utils as ga

# Genetic algorithm
def run_GA(population_size, probe_params):
    population = ga.generate_initial_population(2**population_size, probe_params)
    while len(population) > 1:
        population = ga.calculate_population_fitness(population)
        population = ga.generate_next_population(population)
    population = ga.calculate_population_fitness(population)
    return population

# Hill climbing
def run_Hill_Climbing(population):
    best_probe = population[0]
    population = ga.hill_climb(best_probe)
    population = ga.calculate_population_fitness(population)
    population.sort(reverse=True)
    curr_probe = population[0]
    while curr_probe.beta[0] > best_probe.beta[0]:
        best_probe = curr_probe
        population = ga.hill_climb(best_probe)
        population = ga.calculate_population_fitness(population)
        population.sort(reverse=True)
        curr_probe = population[0]
    sequences = best_probe.sequences
    return sequences

# Clean up output sequences
def handle_output(sequences):
    if len(sequences['sinkC']) < 7:
        sequences['sinkC'] = 'None Required'
        if len(sequences['sink']) < 6:
            sequences['sink'] = 'None Required'
    return sequences

# Run full optimization workflow
def run(population_size, probe_params):
    population = run_GA(population_size, probe_params)
    sequences = run_Hill_Climbing(population)
    final_sequences = handle_output(sequences)
    print(final_sequences)
    print('Done')
    return final_sequences