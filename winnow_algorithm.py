import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# Descarga el archivo de la base de datos
url = 'http://brl.bcm.tmc.edu/rme/gbm.dat'
data = pd.read_csv(url, sep='\t')

# Filtra los genes con al menos 7 muestras mutadas
filtered_data = data.loc[:, (data != 0).sum() > 6]

# Calcula los puntajes de exclusividad entre cada par de genes
exclusivity_scores = {}
for gene1 in filtered_data.columns[1:]:
    for gene2 in filtered_data.columns[1:]:
        if gene1 == gene2:
            continue
        num_samples_where_both_mutated = ((filtered_data[gene1] > 0) & (filtered_data[gene2] > 0)).sum()
        num_samples_where_either_mutated = ((filtered_data[gene1] > 0) | (filtered_data[gene2] > 0)).sum()
        if num_samples_where_either_mutated != 0:
            exclusivity_scores[(gene1, gene2)] = num_samples_where_both_mutated / num_samples_where_either_mutated

# Filtra los puntajes de exclusividad para aquellos que superan un umbral
threshold = 0.8
filtered_exclusivity_scores = {pair: score for pair, score in exclusivity_scores.items() if score > threshold}

# Crea un gráfico de red utilizando NetworkX
G = nx.Graph()
for pair, score in filtered_exclusivity_scores.items():
    gene1, gene2 = pair
    G.add_edge(gene1, gene2, weight=score)

# Visualiza la red de exclusividad
plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_size=500, font_size=10, font_color='black')
labels = nx.get_edge_attributes(G, 'weight')
nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
plt.title('Red de Exclusividad entre Genes')
plt.show()
