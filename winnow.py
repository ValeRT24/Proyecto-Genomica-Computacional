# Aplicar algoritmo de winnow a un conjunto de datos
# Uso: python winnow.py <archivo de datos> <archivo de pesos> <archivo de salida>
# Ejemplo: python winnow.py data.txt weights.txt output.txt

import sys
import numpy as np

# Leer archivo de datos
def read_data(filename):
    data = []
    with open(filename, 'r') as file:
            for line in file:
                data.append(list(map(int, line.strip().split())))
    return data

# Leer archivo de pesos
def read_weights(filename):
    with open(filename, 'r') as file:
        return list(map(float, file.readline().strip().split()))
    
# Escribir archivo de salida
def write_output(filename, data):
    with open(filename, 'w') as file:
        for line in data:
            file.write(str(line) + '\n')

# Algoritmo de winnow
def winnow(data, weights):
    output = []
    for i in range(len(data)):
        x = data[i][:-1]
        y = data[i][-1]
        if np.dot(x, weights) > 0:
            output.append(1)
        else:
            output.append(0)
    return output

# Main
if __name__ == '__main__':
    data = read_data(sys.argv[1])
    weights = read_weights(sys.argv[2])
    output = winnow(data, weights)
    write_output(sys.argv[3], output)
    print('Archivo de salida generado:', sys.argv[3])