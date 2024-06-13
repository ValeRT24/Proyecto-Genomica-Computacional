import matplotlib.pyplot as plt
import numpy as np

# Resultados después de aplicar el algoritmo de winnow
pruebas_50 = {('TCGA-02-0007', 'TCGA-06-0156'): 0.45454545454545453, ('TCGA-02-0007', 'TCGA-06-0219'): 0.5, ('TCGA-02-0016', 'TCGA-02-0048'): 0.4666666666666667, ('TCGA-02-0016', 'TCGA-06-0138'): 0.48936170212765956, ('TCGA-02-0016', 'TCGA-06-0156'): 0.6744186046511628, ('TCGA-02-0016', 'TCGA-06-0219'): 0.5510204081632653, ('TCGA-02-0016', 'TCGA-06-0241'): 0.45454545454545453, ('TCGA-02-0016', 'TCGA-06-0646'): 0.4772727272727273, ('TCGA-02-0016', 'TCGA-08-0347'): 0.475, ('TCGA-02-0023', 'TCGA-08-0348'): 0.5625, ('TCGA-02-0023', 'TCGA-08-0354'): 0.47619047619047616, ('TCGA-02-0023', 'TCGA-08-0359'): 0.47368421052631576, ('TCGA-02-0024', 'TCGA-06-0241'): 0.4772727272727273, ('TCGA-02-0037', 'TCGA-06-0648'): 0.45454545454545453, ('TCGA-02-0048', 'TCGA-02-0016'): 0.4666666666666667, ('TCGA-02-0048', 'TCGA-06-0138'): 0.55, ('TCGA-02-0048', 'TCGA-06-0156'): 0.5609756097560976, ('TCGA-02-0048', 'TCGA-06-0219'): 0.4782608695652174, ('TCGA-02-0048', 'TCGA-06-0241'): 0.6, ('TCGA-02-0048', 'TCGA-06-0646'): 0.6764705882352942, ('TCGA-02-0048', 'TCGA-08-0347'): 0.5, ('TCGA-02-0052', 'TCGA-06-0648'): 0.5555555555555556, ('TCGA-02-0116', 'TCGA-08-0358'): 0.45161290322580644, ('TCGA-02-0116', 'TCGA-08-0375'): 0.4838709677419355, ('TCGA-06-0138', 'TCGA-02-0016'): 0.48936170212765956, ('TCGA-06-0138', 'TCGA-02-0048'): 0.55, ('TCGA-06-0138', 'TCGA-06-0156'): 0.5454545454545454, ('TCGA-06-0138', 'TCGA-06-0219'): 0.46938775510204084, ('TCGA-06-0138', 'TCGA-06-0241'): 0.5, ('TCGA-06-0138', 'TCGA-06-0646'): 0.6052631578947368, ('TCGA-06-0138', 'TCGA-08-0347'): 0.4864864864864865, ('TCGA-06-0150', 'TCGA-06-0211'): 0.5, ('TCGA-06-0150', 'TCGA-08-0358'): 0.5217391304347826, ('TCGA-06-0154', 'TCGA-06-0187'): 0.52, ('TCGA-06-0154', 'TCGA-08-0353'): 0.4642857142857143, ('TCGA-06-0154', 'TCGA-08-0358'): 0.6666666666666666, ('TCGA-06-0154', 'TCGA-08-0375'): 0.5185185185185185, ('TCGA-06-0156', 'TCGA-02-0007'): 0.45454545454545453, ('TCGA-06-0156', 'TCGA-02-0016'): 0.6744186046511628, ('TCGA-06-0156', 'TCGA-02-0048'): 0.5609756097560976, ('TCGA-06-0156', 'TCGA-06-0138'): 0.5454545454545454, ('TCGA-06-0156', 'TCGA-06-0188'): 0.4883720930232558, ('TCGA-06-0156', 'TCGA-06-0219'): 0.574468085106383, ('TCGA-06-0156', 'TCGA-06-0241'): 0.55, ('TCGA-06-0156', 'TCGA-06-0646'): 0.5365853658536586, ('TCGA-06-0156', 'TCGA-08-0347'): 0.5, ('TCGA-06-0158', 'TCGA-08-0357'): 0.5, ('TCGA-06-0159', 'TCGA-08-0244'): 0.47619047619047616, ('TCGA-06-0185', 'TCGA-08-0244'): 0.5185185185185185, ('TCGA-06-0187', 'TCGA-06-0154'): 0.52, ('TCGA-06-0187', 'TCGA-08-0353'): 0.5416666666666666, ('TCGA-06-0187', 'TCGA-08-0355'): 0.47619047619047616, ('TCGA-06-0187', 'TCGA-08-0358'): 0.5652173913043478, ('TCGA-06-0187', 'TCGA-08-0375'): 0.5416666666666666, ('TCGA-06-0188', 'TCGA-06-0156'): 0.4883720930232558, ('TCGA-06-0188', 'TCGA-06-0219'): 0.5454545454545454, ('TCGA-06-0195', 'TCGA-06-0211'): 0.45454545454545453, ('TCGA-06-0211', 'TCGA-06-0150'): 0.5, ('TCGA-06-0211', 'TCGA-06-0195'): 0.45454545454545453, ('TCGA-06-0214', 'TCGA-08-0358'): 0.4642857142857143, ('TCGA-06-0219', 'TCGA-02-0007'): 0.5, ('TCGA-06-0219', 'TCGA-02-0016'): 0.5510204081632653, ('TCGA-06-0219', 'TCGA-02-0048'): 0.4782608695652174, ('TCGA-06-0219', 'TCGA-06-0138'): 0.46938775510204084, ('TCGA-06-0219', 'TCGA-06-0156'): 0.574468085106383, ('TCGA-06-0219', 'TCGA-06-0188'): 0.5454545454545454, ('TCGA-06-0219', 'TCGA-06-0646'): 0.4888888888888889, ('TCGA-06-0219', 'TCGA-08-0347'): 0.4523809523809524, ('TCGA-06-0241', 'TCGA-02-0016'): 0.45454545454545453, ('TCGA-06-0241', 'TCGA-02-0024'): 0.4772727272727273, ('TCGA-06-0241', 'TCGA-02-0048'): 0.6, ('TCGA-06-0241', 'TCGA-06-0138'): 0.5, ('TCGA-06-0241', 'TCGA-06-0156'): 0.55, ('TCGA-06-0241', 'TCGA-06-0646'): 0.5714285714285714, ('TCGA-06-0241', 'TCGA-08-0347'): 0.53125, ('TCGA-06-0646', 'TCGA-02-0016'): 0.4772727272727273, ('TCGA-06-0646', 'TCGA-02-0048'): 0.6764705882352942, ('TCGA-06-0646', 'TCGA-06-0138'): 0.6052631578947368, ('TCGA-06-0646', 'TCGA-06-0156'): 0.5365853658536586, ('TCGA-06-0646', 'TCGA-06-0219'): 0.4888888888888889, ('TCGA-06-0646', 'TCGA-06-0241'): 0.5714285714285714, ('TCGA-06-0646', 'TCGA-08-0347'): 0.5151515151515151, ('TCGA-06-0648', 'TCGA-02-0037'): 0.45454545454545453, ('TCGA-06-0648', 'TCGA-02-0052'): 0.5555555555555556, ('TCGA-08-0244', 'TCGA-06-0159'): 0.47619047619047616, ('TCGA-08-0244', 'TCGA-06-0185'): 0.5185185185185185, ('TCGA-08-0244', 'TCGA-08-0354'): 0.5416666666666666, ('TCGA-08-0347', 'TCGA-02-0016'): 0.475, ('TCGA-08-0347', 'TCGA-02-0048'): 0.5, ('TCGA-08-0347', 'TCGA-06-0138'): 0.4864864864864865, ('TCGA-08-0347', 'TCGA-06-0156'): 0.5, ('TCGA-08-0347', 'TCGA-06-0219'): 0.4523809523809524, ('TCGA-08-0347', 'TCGA-06-0241'): 0.53125, ('TCGA-08-0347', 'TCGA-06-0646'): 0.5151515151515151, ('TCGA-08-0348', 'TCGA-02-0023'): 0.5625, ('TCGA-08-0348', 'TCGA-08-0354'): 0.45454545454545453, ('TCGA-08-0348', 'TCGA-08-0359'): 0.8125, ('TCGA-08-0353', 'TCGA-06-0154'): 0.4642857142857143, ('TCGA-08-0353', 'TCGA-06-0187'): 0.5416666666666666, ('TCGA-08-0353', 'TCGA-08-0358'): 0.5, ('TCGA-08-0353', 'TCGA-08-0375'): 0.48148148148148145, ('TCGA-08-0354', 'TCGA-02-0023'): 0.47619047619047616, ('TCGA-08-0354', 'TCGA-08-0244'): 0.5416666666666666, ('TCGA-08-0354', 'TCGA-08-0348'): 0.45454545454545453, ('TCGA-08-0355', 'TCGA-06-0187'): 0.47619047619047616, ('TCGA-08-0355', 'TCGA-08-0356'): 0.5555555555555556, ('TCGA-08-0356', 'TCGA-08-0355'): 0.5555555555555556, ('TCGA-08-0357', 'TCGA-06-0158'): 0.5, ('TCGA-08-0358', 'TCGA-02-0116'): 0.45161290322580644, ('TCGA-08-0358', 'TCGA-06-0150'): 0.5217391304347826, ('TCGA-08-0358', 'TCGA-06-0154'): 0.6666666666666666, ('TCGA-08-0358', 'TCGA-06-0187'): 0.5652173913043478, ('TCGA-08-0358', 'TCGA-06-0214'): 0.4642857142857143, ('TCGA-08-0358', 'TCGA-08-0353'): 0.5, ('TCGA-08-0358', 'TCGA-08-0375'): 0.56, ('TCGA-08-0359', 'TCGA-02-0023'): 0.47368421052631576, ('TCGA-08-0359', 'TCGA-08-0348'): 0.8125, ('TCGA-08-0375', 'TCGA-02-0116'): 0.4838709677419355, ('TCGA-08-0375', 'TCGA-06-0154'): 0.5185185185185185, ('TCGA-08-0375', 'TCGA-06-0187'): 0.5416666666666666, ('TCGA-08-0375', 'TCGA-08-0353'): 0.48148148148148145, ('TCGA-08-0375', 'TCGA-08-0358'): 0.56}
pruebas_60 = {('TCGA-02-0016', 'TCGA-06-0156'): 0.6744186046511628, ('TCGA-02-0048', 'TCGA-06-0646'): 0.6764705882352942, ('TCGA-06-0154', 'TCGA-08-0358'): 0.6666666666666666, ('TCGA-06-0156', 'TCGA-02-0016'): 0.6744186046511628, ('TCGA-06-0646', 'TCGA-02-0048'): 0.6764705882352942, ('TCGA-08-0348', 'TCGA-08-0359'): 0.8125, ('TCGA-08-0358', 'TCGA-06-0154'): 0.6666666666666666, ('TCGA-08-0359', 'TCGA-08-0348'): 0.8125}

# Para eliminar pares invertidos
def eliminar_pares_invertidos(diccionario):
    nuevo_diccionario = {}
    for key, value in diccionario.items():
        if (key[1], key[0]) in nuevo_diccionario:
            continue
        nuevo_diccionario[key] = value
    return nuevo_diccionario

pruebas_5 = eliminar_pares_invertidos(pruebas_50)
pruebas_6 = eliminar_pares_invertidos(pruebas_60)

# Función que nos permite calcula la significancia RME con el teorema de Bayes.
def valor_significancia_RME(diccionario, P_xij_1=0.5, P_xij_0=0.5):
    valores = {}
    for (gen_a, gen_b), P_xij_ab in diccionario.items():
        # probabilidades condicionadas P(x_ij = 1 | a_ij)
        P_xij_a = P_xij_ab  
        # probabilidades condicionadas P(x_ij = 1 | b_ij)
        P_xij_b = P_xij_ab  

        # Usando la fórmula del teorema de Bayes modificada:
        numerador = P_xij_a * P_xij_b * P_xij_1
        denominador = ((1 - P_xij_a) * (1 - P_xij_b) * P_xij_0 + P_xij_a * P_xij_b * P_xij_1)
        rme = numerador / denominador

        # Almacenamos el resultado
        valores[(gen_a, gen_b)] = rme
    
    return valores

# Función para calcular d_prime usando la matriz de aberración
def codificar_matriz_aberracion(matriz, probabilidades_nulas, probabilidades_rme):
    d = 0
    num_filas, num_columnas = matriz.shape
    
    for i in range(num_filas):
        for j in range(num_columnas):
            par_genes = (f'TCGA-02-00{str(i+1).zfill(2)}', f'TCGA-06-01{str(j+1).zfill(2)}')
            if par_genes in probabilidades_nulas and par_genes in probabilidades_rme:
                p_null = probabilidades_nulas[par_genes]
                p_rme = probabilidades_rme[par_genes]
                if matriz[i, j] == 1:
                    d += (-np.log2(p_null)) + np.log2(p_rme)
                else:
                    d += (-np.log2(p_null)) + np.log2(1 - p_rme)
    return d

# Para aplicar el algoritmo nulo y calcular d
def algoritmo_nulo(diccionario, probabilidades_nulas, probabilidades_rme):
    d = 0
    for (genes_nulos, genes_rme) in zip(probabilidades_nulas, probabilidades_rme):
        p_null = probabilidades_nulas[genes_nulos]
        p_rme = probabilidades_rme[genes_rme]
        if p_null >= 0.5:
            d += (-np.log2(p_null)) + np.log2(p_rme)
        else:
            d += (-np.log2(p_null)) + np.log2(1 - p_rme)
    
    return d


# Calcular la significancia final
def calcular_significancia(d_valor, m, k):
    d = d_valor - (m * np.log2(k) - k * np.log2(m) - m * np.log2(k))
    significancia = 2 ** (-d)
    return significancia

# Función para eliminar los genes no mutados
def eliminar_genes_no_mutados(diccionario1, diccionario2):
    # Crear una lista de claves para evitar modificar el diccionario mientras lo recorremos
    claves_a_eliminar = []
    
    for (gen1_a, gen1_b), valor1 in diccionario1.items():
        for (gen2_a, gen2_b), valor2 in diccionario2.items():
            # Condición para determinar si un gen debe ser eliminado
            if valor1 > valor2:
                claves_a_eliminar.append((gen2_a, gen2_b))
    
    # Eliminar las claves marcadas
    for clave in claves_a_eliminar:
        diccionario2.pop(clave, None)
    
    return diccionario2

# Ejecutar el algoritmo a través del resultado del algoritmo de winnow
def obtener_significancia(diccionario):
    # Ejecutar el algoritmo
    valores_significancia = valor_significancia_RME(diccionario)
    # Algoritmo nulo 
    d = algoritmo_nulo(diccionario, diccionario, valores_significancia)

    # Matriz
    m = len(diccionario)  
    k = len(diccionario)

    # Calcular la significancia en general
    significancia_final = calcular_significancia(d, m, k)   

    print(significancia_final)

    # Elminar los genes no mutados
    valores_significancia_final = eliminar_genes_no_mutados(diccionario, valores_significancia)
    
    for genes, significancia in valores_significancia_final.items():
        print(f"Porcentaje de Significancia para {genes}: {significancia:.4f}")

    return valores_significancia

diccionario1 = pruebas_5
diccionario2 = valor_significancia_RME(pruebas_5)

obtener_significancia(pruebas_5)

# Comparar los valores y calcular diferencias
diferencias = {}
for key in diccionario1:
    if key in pruebas_5:
        diferencias[key] = diccionario2[key] - pruebas_5[key]

# Crear listas para los nombres de pares y las diferencias
pares = list(diferencias.keys())
valores_diferencia = list(diferencias.values())

# Convertir los nombres de pares en strings para usarlos en las etiquetas del gráfico
etiquetas = [f"{gen_a}\n{gen_b}" for gen_a, gen_b in pares]

# Configurar el gráfico
fig, ax = plt.subplots(figsize=(10, 6))

# Crear una barra para cada par de genes
bar_width = 0.35
indices = np.arange(len(pares))
barras = ax.bar(indices, valores_diferencia, bar_width, color='r')

# Añadir etiquetas y títulos
ax.set_xlabel('Pares de genes')
ax.set_ylabel('Diferencia de valores')
ax.set_title('Porcentaje de aumento en mutación de los genes')
ax.set_xticks(indices)
ax.set_xticklabels(etiquetas, rotation=45, ha='right')
ax.axhline(0, color='grey', linewidth=0.8)

# Añadir etiquetas de valor en las barras
for barra in barras:
    altura = barra.get_height()
    ax.text(barra.get_x() + barra.get_width() / 2.0, altura, f'{altura:.2f}', ha='center', va='bottom')

# Mostrar el gráfico
plt.tight_layout()
plt.show()
