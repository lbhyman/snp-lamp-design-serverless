
def zeros_matrix(n_rows, n_cols):
    output = []
    while len(output) < n_rows:
        output.append([0.0]*n_cols)
    return output

def transpose(mat):
    n_rows, n_cols = len(mat), len(mat[0])
    output = zeros_matrix(n_cols, n_rows)
    for i in range(n_rows):
        for j in range(n_cols):
            output[j][i] = mat[i][j]
    return output

def dot(mat_1, mat_2):
    output = zeros_matrix(len(mat_1), len(mat_2[0]))
    for i in range(len(mat_1)):
        for j in range(len(mat_2[0])):
            total = 0
            for ii in range(len(mat_1[0])):
                total += mat_1[i][ii] * mat_2[ii][j]
            output[i][j] = total
    return output