from chalice import Chalice
from probe import Probe
import ga_utils as ga
import json

app = Chalice(app_name='snp-lamp-design')

# Generate initial GA population
@app.route('/generate_initial_population/{params}')
def generate_initial_population(params):
    params = json.loads(params)
    pop_size = int(params['popSize'])
    probe_params = params['probeParams']
    return json.dumps(ga.generate_initial_population(pop_size,probe_params))

# Generate next GA population
@app.route('/generate_next_population/{population}')
def generate_next_population(population):
    return json.dumps(ga.generate_next_population(json.loads(population)))

# Generate next hill-climbing steps
@app.route('/hill_climbing_options/{probe_list}')
def hill_climb(probe_list):
    new_list = ga.hill_climb(json.loads(probe_list))
    return json.dumps(new_list)

# Calculate probe fitness
@app.route('/calculate_fitness/{probe}')
def calculate_fitness(probe):
    curr_probe = Probe(json.loads(probe))
    curr_probe.calc_beta()
    return json.dumps(curr_probe.__dict__)
