import chalicelib.LA_utils as la
from Levenshtein import distance
import boto3
from boto3.dynamodb.conditions import Key

BASE_ENCODINGS = {
    'A': [1.0,0.0,0.0,0.0],
    'T': [0.0,1.0,0.0,0.0],
    'G': [0.0,0.0,1.0,0.0],
    'C': [0.0,0.0,0.0,1.0]
}

def get_db():
    dynamodb = boto3.resource("dynamodb")
    table = dynamodb.Table('snp-lamp-design-table-lsh')
    return table

class LSH():
    
    def __init__(self, num_tables, num_bits, input_dim, projection_filename):
        self.input_dim = input_dim*4
        self.num_bits = num_bits
        self.num_tables = num_tables
        self.projections = self.load_projections(projection_filename)
        self.table = get_db()

    def get(self, item):
        results = []
        hash_vals = self.hash(item)
        '''response = self.table.query(
            IndexName='LSH',
            KeyConditionExpression=Key('hash_loc').eq(hash_vals[0]) & 
            (
                Key('hash_1').eq(hash_vals[1]) | 
                Key('hash_2').eq(hash_vals[2]) | 
                Key('hash_3').eq(hash_vals[3]) | 
                Key('hash_4').eq(hash_vals[4])
            )
        )'''
        response = self.table.query(
            IndexName='hash_loc-index',
            KeyConditionExpression=Key('hash_loc').eq(hash_vals[0])# & Key('hash_1').eq(hash_vals[1])
        )
        for r in response['Items']:
            curr_item = r
            if ((curr_item['hash_1'] == hash_vals[1]) |
                (curr_item['hash_2'] == hash_vals[2]) |
                (curr_item['hash_3'] == hash_vals[3]) |
                (curr_item['hash_4'] == hash_vals[4])):
                curr_item['distance'] = distance(item['WT'], r['seq'])
                results.append(curr_item)
        sorted_results = sorted(results, key= lambda x: x['distance'])
        return sorted_results
    
    def put(self, id, item):
        hash_vals = self.hash(item)
        item['hash_loc'] = hash_vals[0]
        item['hash_1'] = hash_vals[1]
        item['hash_2'] = hash_vals[2]
        item['hash_3'] = hash_vals[3]
        item['hash_4'] = hash_vals[4]
        table.put_item(Item=item)
    
    def load_projections(self, filename):
        projections, curr_projection = [], []
        infile = open(filename, 'r')
        for line in infile.readlines():
            if '>' in line:
                projections.append(curr_projection)
                curr_projection = []
            else:
                line = line.strip().split(',')
                curr_projection.append([float(val) for val in line])
        infile.close()
        return projections
            
    def encode_len(self, seq):
        length = 0
        if len(seq) <= 3:
            length = len(seq)
        else:
            # round to nearest three bases
            length = 3 * int(float(len(seq))/3.0) 
        binary = '{0:b}'.format(length)
        return ('0' * (8 - len(binary))) + binary

    def encode_seq(self, seq, projection):
        output = []
        for i in range(len(seq)):
            output = output + BASE_ENCODINGS[seq[i]]
        output = output + [0.0]*(self.input_dim - len(output))
        output = la.dot([output], la.transpose(projection))[0]
        final_output = ''
        for i in range(len(output)):
            if output[i] > 0:
                final_output = final_output + '1'
            else:
                final_output = final_output + '0'
        return final_output
            
    def hash(self, item):
        WT, SNP_index = item['WT'], item['SNP_index']
        hash_vals = [self.encode_len(WT[:SNP_index]) + self.encode_len(WT[SNP_index+1:])]
        for projection in self.projections:
            hash_vals += [self.encode_seq(WT, projection)]
        return hash_vals
