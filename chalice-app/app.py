from chalice import Chalice, Response, CORSConfig
from chalicelib.probe import Probe
import chalicelib.ga_utils as ga
from chalicelib import optimizer
import json
import boto3
from boto3.dynamodb.conditions import Key

app = Chalice(app_name='snp-lamp-design')

cors_config = CORSConfig(
    allow_origin='http://192.168.0.158:3000',
    allow_credentials=True
)

def get_app_db():
    dynamodb = boto3.resource("dynamodb")
    table = dynamodb.Table('snp-lamp-design-table')
    return table

def generate_jobID(params):
    probeParams = params['probeParams']
    jobID = '/'.join([probeParams['WT'], probeParams['SNP'], str(params['popSize']), str(probeParams['params']['temperature']), str(probeParams['params']['magnesium']), str(probeParams['params']['sodium'])])
    return jobID

@app.route('/', methods=['GET', 'POST'], cors=True)
def index():
    request = app.current_request
    return json.dumps({'hello': 'world'})

# Start Optimizer
@app.route('/start_optimizer', methods=['POST'], cors=True)
def start_optimizer():
    request = app.current_request
    if request.method == 'POST':
        params = json.loads(request._body)
        pop_size = int(params['popSize'])
        probe_params = params['probeParams']
        jobID = generate_jobID(params)
        table = get_app_db()
        # Check if job already submitted
        response = table.get_item(Key={'id': jobID, 'SNP': probe_params['SNP']})
        if response.has_key('Item'):
            item = response['Item']
            return json.dumps(item)
        # Start job
        else:
            table.put_item(Item={
                'id': jobID,
                'SNP': probe_params['SNP'],
                'WT': probe_params['WT'],
                'probeF': '',
                'probeQ': '',
                'sink': '',
                'sinkC': '',
                'running': 'true'
            })
            output_probe = optimizer.run(pop_size, probe_params)
            output_probe['id'] = jobID
            output_probe['running'] = 'false'
            table.put_item(Item=output_probe)
            return json.dumps(output_probe)

# Get Output
@app.route('/get_output', methods=['POST'], cors=True)
def get_output():
    request = app.current_request
    if request.method == 'POST':
        params = json.loads(request._body)
        jobID = generate_jobID(params)
        table = get_app_db()
        response = table.get_item(Key={'id': jobID, 'SNP': probe_params['SNP']})
        item = response['Item']
        return json.dumps(item)