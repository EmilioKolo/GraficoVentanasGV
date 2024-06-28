
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

###########################################################################################
#########################################FUNCIONES#########################################
###########################################################################################

def abrirExcel(nombre,L_head,extension='xlsx'):
    #Parseo del nombre del archivo
    file_name = os.path.dirname(os.path.abspath(__file__)) + '\\' + str(nombre) + '.' + str(extension);

    #Abro el archivo con Pandas
    df = pd.read_excel(io=file_name);

    L = [];

    #Hago una lista de cada columna en L_heads
    for i in L_head:
        L.append(list(df[i]));

    return L

###########################################################################################
#########################################EJECUCION#########################################
###########################################################################################

#Nombre del archivo
nom = 'Fig3';
#Nombres de columnas
L_head = ['mRNA','Pathogenic','Healthy','Protein'];

#Creacion de lista de columnas
L_R = abrirExcel(nom,L_head);

#Asignacion de columnas a sus variables
mRNA, Path, Hthy, Prot = L_R[0], L_R[1], L_R[2], L_R[3];

###########################################################################################
#########################################GRAFICOS##########################################
###########################################################################################

#Defino el eje x
x = Prot;
xlim = [0,325]

#Listas de tamaños
T0 = (12,4,3,8);
T1 = (24,8,4,10);
T = T0;

#Ancho de linea
lw = 1;

#Tamaño total de la figura
fsize = (T[0],T[1]);

plt.figure(num=None, figsize=fsize, dpi=80, facecolor='w', edgecolor='k')

#Armo el plot doble
plt.subplot(T[2],1,(1,T[2]-1));

#Grafico ambos sets de datos
#plt.plot(x,Hthy,linewidth=lw,color='g',label='Frequency>1/10,000');
plt.plot(x,Path,linewidth=lw,color='r',label='Pathogenic');

#Defino el diseño de los ejes
plt.xlim(xlim);
plt.xticks(range(0,xlim[1]+25,25));
plt.xlabel('Position in NKX2-5 protein');
plt.ylim([-0.8,7.5]);
plt.yticks(range(0,8,1));
plt.ylabel('Number of pathogenic GVs')

#Muestro la leyenda
#plt.legend(loc='upper right')

##Defino el titulo
#plt.title('Figure 3');

#Trabajo sobre la segunda parte
plt.subplot(T[2],1,T[2]);

#Ancho de las regiones/dominios
dl = T[3];

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

'''
color='blue'        # specify color by name
color='g'           # short color code (rgbcmyk)
color='0.75'        # Grayscale between 0 and 1
color='#FFDD44'     # Hex code (RRGGBB from 00 to FF)
color=(1.0,0.2,0.3) # RGB tuple, values 0 to 1
color='chartreuse'  # all HTML color names supported
'''

#Escondo los ejes
plt.axis('off');

#Muestro o guardo el grafico
plt.savefig('Fig3.png', dpi=100);
print('Figura guardada.');
#plt.show();
'''
# (1) save the image in memory in PNG format
png1 = BytesIO()
plt.savefig(png1, format='png')

# (2) load this image into PIL
png2 = Image.open(png1)

# (3) save as TIFF
png2.save('3dPlot.tiff')
png1.close()
'''
'''
#Muestro o guardo la figura
b = str(input('0-Salir\n1-Guardar la figura\nCualquier otra opcion hace display.\n'));
if b == '1':
    plt.savefig('Fig3.png', dpi=100);
    print('Figura guardada.');
elif b != '0':
    plt.show();
else:
    plt.clf();
    exit()
'''
