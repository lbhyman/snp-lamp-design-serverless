from chalice import Chalice
from chalice import CORSConfig
from probe import Probe
import ga_utils as ga
import optimizer
import json
import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=FutureWarning)

app = Chalice(app_name='snp-lamp-design')

cors_config = CORSConfig(
    allow_origin='http://localhost:3000',
    allow_credentials=True
)

# Run Optimizer
@app.route('/start_optimizer', methods=['POST'], cors=cors_config)
def start_optimizer():
    request = app.current_request
    if request.method == 'POST':
        params = json.loads(request._body)
        pop_size = int(params['popSize'])
        probe_params = params['probeParams']
        output_probe = optimizer.run(pop_size, probe_params)
        return json.dumps(output_probe)
