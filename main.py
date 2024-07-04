
# Generales
import os 
import time 
import copy 
import math 
import sys 
import pandas as pd 
import numpy as np 

# Para display
from random import shuffle 

# Graficos
import matplotlib.pyplot as plt 
import seaborn as sns 
import pandas as pd


'''
Funciones para generar un grafico de cantidad de variantes en una ventana de nucleotidos
'''

#################################### VARIABLES ####################################

### Variables main()

# Path main depende de posicion actual
path_git_main = os.path.dirname(os.path.abspath(__file__)); 
path_output_dump_main = 'D:\\ArchivosDoctorado\\Output_dump\\'; 
path_in_main = path_git_main + ''; ### COMPLETAR
path_out_main = path_output_dump_main + ''; ### COMPLETAR

nom_archivo_variantes = 'tabla_variantes'; 
path_archivo_variantes = path_git_main; 

# Columnas seleccionadas, las primeras 2 son identificadores, 
# la 3ra es para seleccionar pathogenic y healthy, la 4ta es para seleccionar missense
# la 5ta es para seleccionar frecuencias
nom_columnas = ["cDNA", "Protein", "Classification", "Protein effect", "Frequency per 1000 ind"]; 

nom_out_main = 'FiguraVentanas'; 


#################################### FUNCIONES ####################################


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    pipeline_grafico(nom_archivo_variantes, nom_out_main, nom_columnas, path_gvs=path_archivo_variantes, path_out=path_archivo_variantes); 

    return ''


def pipeline_grafico(nom_gvs, nom_out, nom_cols, path_gvs='', path_out=''): 
    '''Genera el grafico de GVs por ventana desde el archivo de variantes y las columnas correspondientes'''

    # Abro el archivo y extraigo la columna
    m_gvs, l_head = abrir_csv(nom_gvs, path_arch=path_gvs, con_headers=True, devolver_header=True); 

    # Uso la funcion seleccionar_columnas_filtro para seleccionar los datos correctos y transformarlos
    m_gvs_p, m_gvs_b = seleccionar_columnas_filtro(m_gvs, l_head, nom_cols); 

    # Uso la funcion generar_grafico para hacer graficos con m_gvs_p y m_gvs_b
    generar_grafico(m_gvs_p, m_gvs_b, nom_out, path_out=path_out); 
    
    return m_gvs_p, m_gvs_b


### Funciones principales del pipeline

def contar_para_grafico(m_gvs_p, m_gvs_b, rango=(1,975), ventana=9):
    '''Cuenta la cantidad de variantes en ventanas dentro del rango para m_gvs_b y m_gvs_p'''

    # Inicializo las matrices que se devuelven
    x = []; 
    col_p = []; 
    col_b = []; 
    # Recorro el rango con la ventana
    for i in range(rango[0], rango[1]-ventana+2):
        # Defino posiciones iniciales, finales y el centro
        pos_ini = i; 
        pos_end = i + ventana - 1; 
        medio = (pos_ini+pos_end)/2; 
        # Agrego medio a x
        x.append(float(medio)/3); 
        # Uso contar_hits() para ver cuantos hits hay por posicion
        col_p.append(contar_hits(m_gvs_p, pos_ini, pos_end, col_id=0)); 
        col_b.append(contar_hits(m_gvs_b, pos_ini, pos_end, col_id=0)); 
    return x, col_p, col_b


def generar_grafico(m_gvs_p, m_gvs_b, nom_out, path_out=''):
    '''Funcion que hace un grafico en base a los datos ya procesados'''

    xlim = [0,325]; 
    # Paso los datos a listas por columnas
    x, col_n_path, col_n_benign = contar_para_grafico(m_gvs_p, m_gvs_b); 
    # Listas de tamaños
    t0 = (12,4,3,8); 
    t1 = (24,8,4,10); 
    t = t0; 
    # Ancho de linea
    lw = 1; 
    # Tamaño total de la figura
    fsize = (t[0],t[1]); 
    plt.figure(num=None, figsize=fsize, dpi=80, facecolor='w', edgecolor='k'); 
    # Armo el plot doble
    plt.subplot(t[2],1,(1,t[2]-1)); 
    #Grafico ambos sets de datos
    plt.plot(x,col_n_benign,linewidth=lw,color='g',label='Frecuencia mayor a 1/10.000'); 
    plt.plot(x,col_n_path,linewidth=lw,color='r',label='Variantes patogénicas'); 

    #Defino el diseño de los ejes
    plt.xlim(xlim); 
    plt.xticks(range(0,xlim[1]+25,25)); 
    plt.xlabel('Posición en la proteína NKX2-5'); 
    plt.ylim([-0.8,7.5]); 
    plt.yticks(range(0,8,1)); 
    plt.ylabel('Número de variantes'); 

    #Muestro la leyenda
    plt.legend(loc='upper right')
    
    #Trabajo sobre la segunda parte
    plt.subplot(t[2],1,t[2]); 

    #Ancho de las regiones/dominios
    dl = t[3]; 

    #Posicion de las regiones/dominios
    plt.ylim([-100,100]); 
    plt.yticks(range(-100,125,25)); 
    h = [0, 0]; 
    dif = dl*2; 

    #mRNA
    TN = [10*3-2+dl/2,21*3-dl/2]; 
    HD = [138*3-2+dl/2,197*3-dl/2]; 
    H1 = [146*3-2+dl/2,158*3-dl/2]; 
    H2 = [165*3-2+dl/2,176*3-dl/2]; 
    H3 = [179*3-2+dl/2,194*3-dl/2]; 
    NK2SD = [212*3-2+dl/2,234*3-dl/2]; 
    YRR = [237*3-2+dl/2,274*3-dl/2]; 
    NKbox = [291*3-2+dl/2,304*3-dl/2]; 
    GIRAW = [320*3-2+dl/2,324*3-dl/2]; 
    #Prot
    TN = [10+dl/6,21-dl/6]; 
    HD = [138+dl/6,197-dl/6]; 
    NK2SD = [212+dl/6,234-dl/6]; 
    YRR = [237+dl/6,274-dl/6]; 
    NKbox = [291+dl/6,304-dl/6]; 
    GIRAW = [320+dl/6,324-dl/6]; 
    H1 = [146+dl/6,158-dl/6]; 
    H2 = [165+dl/6,176-dl/6]; 
    H3 = [179+dl/6,194-dl/6]; 

    #Plot de las regiones/dominios
    plt.plot(TN,h,color='#C8C800',linewidth=dl,fillstyle='none'); #Tinman
    plt.plot(HD,[-dif,-dif],color='#640000',linewidth=dl,fillstyle='none'); #Homeodominio
    plt.plot(H1,[dif,dif],color='r',linewidth=dl,fillstyle='none'); #Helice alfa 1
    plt.plot(H2,[dif,dif],color='r',linewidth=dl,fillstyle='none'); #Helice alfa 2
    plt.plot(H3,[dif,dif],color='r',linewidth=dl,fillstyle='none'); #Helice alfa 3
    plt.plot(NK2SD,h,color='#7D007D',linewidth=dl,fillstyle='none'); #NK2-SD
    plt.plot(YRR,h,color='#00007D',linewidth=dl,fillstyle='none'); #Tyr-Rich region
    plt.plot(NKbox,h,color='#00A2E8',linewidth=dl,fillstyle='none'); #NKX2-5 box
    plt.plot(GIRAW,h,color='#22B14C',linewidth=dl,fillstyle='none'); #GIRAW motif
    plt.xlim(xlim); 

    #Escondo los ejes
    plt.axis('off'); 

    #Muestro o guardo el grafico
    plt.savefig(os.path.join(path_out, nom_out+'.png'), dpi=100); 
    print('Figura guardada.'); 


def seleccionar_columnas_filtro(m_gvs, l_head, nom_cols): 
    '''Selecciona columnas en nom_cols y filtra y procesa las cosas'''

    # Inicializo las matrices que se devuelven
    m_out_pathogenic = []; 
    m_out_benign = []; 
    # Defino ids de nom_cols
    l_ids = []; 
    for k in range(len(nom_cols)): 
        # Agrego el id correspondiente a l_ids
        l_ids.append(l_head.index(nom_cols[k])); 
    # Recorro m_gvs
    for i in range(len(m_gvs)):
        curr_gv = m_gvs[i]; 
        # Inicializo las listas que se agregan
        l_out_p = []; 
        l_out_b = []; 
        # Selecciono missense
        if curr_gv[l_ids[3]].lower()=='missense': 
            # Transformo identificadores en posiciones
            pos_adnc = procesar_adnc(curr_gv[l_ids[0]]); 
            pos_prot = procesar_prot(curr_gv[l_ids[1]]); 
            # Selecciono pathogenic
            if curr_gv[l_ids[2]].lower()=='pathogenic': 
                # Cargo posiciones en l_out_p
                l_out_p.append(pos_adnc); 
                l_out_p.append(pos_prot); 
                l_out_p.append('pathogenic'); 
                # Agrego l_out_p a m_out_pathogenic
                m_out_pathogenic.append(l_out_p[:]); 
            # Selecciono benignas y frecuencias altas
            elif (curr_gv[l_ids[2]].lower() in ['likely benign', 'benign']) or (curr_gv[l_ids[4]]!='' and float(curr_gv[l_ids[4]])>=0.1): 
                # Cargo posiciones en l_out_b
                l_out_b.append(pos_adnc); 
                l_out_b.append(pos_prot); 
                l_out_b.append('benign_freq_alta'); 
                # Agrego l_out_b a m_out_benign
                m_out_benign.append(l_out_b[:]); 
    return m_out_pathogenic, m_out_benign


### Funciones simples de procesamiento


def contar_hits(m_gvs, pos_ini, pos_end, col_id=0):
    '''Cuenta la cantidad de elementos en m_gvs que tienen valor de col_id entre pos_ini y pos_end'''

    # Inicializo el contador que se devuelve
    cont_ret = 0; 
    # Recorro m_gvs
    for i in range(len(m_gvs)):
        curr_gv = m_gvs[i]; 
        # Veo si curr_gv[col_id] esta entre pos_ini y pos_end
        if curr_gv[col_id]>=pos_ini and curr_gv[col_id]<=pos_end:
            cont_ret = cont_ret + 1; 
    return cont_ret


def procesar_adnc(id_adnc):
    '''Recibe un identificador de ADNc en formato HGVS y devuelve solo la posicion'''

    # Inicializo el texto que se devuelve
    str_out = ''; 
    # Borro el c. del principio y el cambio N>N del final
    str_out = id_adnc[2:-3]; 
    return int(str_out)


def procesar_prot(id_prot):
    '''Recibe un identificador de proteina en formato HGVS y devuelve solo la posicion'''

    # Inicializo el texto que se devuelve
    str_out = ''; 
    # Veo si id_prot tiene parentesis
    if "(" in id_prot:
        # Borro el principio y el final
        str_out = id_prot[6:-4]; 
    else:
        # Borro el principio y el final
        str_out = id_prot[5:-3]; 
    return int(str_out)


### Funciones para abrir y guardar archivos

def abrir_csv(nom_arch, path_arch='', sep=';', ext='.csv', con_headers=True, devolver_header=False):
    '''Abre archivos .csv y devuelve una lista de listas con las filas.
    con_headers define si hay header (se ignora por defecto).
    devolver_header permite devolver el header como segundo output'''

    # Inicializo la matriz que se devuelve
    M_out = []; 
    # Inicializo booleano y lista de headers
    header = con_headers; 
    l_head = []; 
    # Defino la direccion del archivo con nom_arch y path_arch
    if path_arch=='':
        filepath = nom_arch + ext; 
    else:
        filepath = os.path.join(path_arch, nom_arch + ext); 
    # Abro el archivo filepath
    with open(filepath, 'r') as f_out:
        print('Archivo ' + nom_arch + ' abierto.')
        # Recorro cada fila
        for curr_line in f_out.readlines():
            # Veo si tengo que pasar el header
            if header:
                # Transformo la linea en lista
                l_head = curr_line.rstrip().split(sep=sep); 
                # Paso header a False para no ver la proxima linea
                header = False; 
            else:
                # Transformo la linea en lista
                l_line = curr_line.rstrip().split(sep=sep); 
                # Cargo l_line en M_out
                M_out.append(l_line[:]); 
    # Si devolver_header es True, se devuelve M_out y l_head
    if devolver_header:
        return M_out, l_head
    return M_out



#################################### RESULTADOS ###################################

output_dump = []; 

if __name__=='__main__':
    output_dump.append(_main()); 

