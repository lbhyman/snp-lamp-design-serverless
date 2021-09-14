from chalice import Chalice, Response, CORSConfig
from chalicelib.probe import Probe
from chalicelib.optimizer import Optimizer
import json

app = Chalice(app_name='snp-lamp-design')

cors_config = CORSConfig(
    allow_origin='http://192.168.0.158:3000',
    allow_credentials=True
)

@app.route('/', methods=['GET', 'POST'], cors=True)
def index():
    request = app.current_request
    return json.dumps({'hello': 'world'})

# Start Optimizer
@app.route('/start_optimizer', methods=['POST'], cors=True)
def start_optimizer():
    request = app.current_request
    if request.method == 'POST':
        body = json.loads(request._body)
        item = {
            'WT': body['WT'],
            'SNP': body['SNP']
        }
        opt = Optimizer()
        sequences = opt.optimize(item)
        if sequences == None:
            return json.dumps({})
        return json.dumps(sequences)