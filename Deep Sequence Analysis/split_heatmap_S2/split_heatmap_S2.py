from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np
from matplotlib.pyplot import figure

# This code only works for 21 x 21 heatmap; must change M and N to modify rows and columns.
# Python file run in JupyterLab

unsorted = 'freq_unsorted.csv'
design = 'freq_design.csv'
save = 'lib_design.png'

M, N = 21, 21  # columns, rows

df_unsorted=pd.read_csv(unsorted, sep=',',header=None)
df_design=pd.read_csv(design, sep=',',header=None)
arr_unsorted = df_unsorted.to_numpy()
arr_design = df_design.to_numpy()
valuesNW = arr_design
valuesSE = arr_unsorted
values = [valuesNW, valuesSE]

xv, yv = np.meshgrid(np.arange(-0.5, M), np.arange(-0.5, N))  # vertices of the little squares
x = xv.ravel()
y = yv.ravel()

trianglesNW = [(i + j*(M + 1), i+1 + j*(M + 1), i + (j + 1)*(M + 1)) for j in range(N) for i in range(M)]
trianglesSE = [(i+1 + j*(M + 1),i+1 + (j + 1)*(M + 1), i + (j + 1) * (M + 1)) for j in range(N) for i in range(M)]
triangul = [Triangulation(x, y, triangles) for triangles in [trianglesNW, trianglesSE]]

plt.rcParams["font.family"] = "Verdana"
plt.rcParams['font.size'] = 11.5

fig, ax = plt.subplots(figsize=(8,8), dpi=150)
imgs = [ax.tripcolor(t, val.ravel(), cmap=cmap, vmin=0, vmax=0.2, ec='black', linewidth=0.01) 
        for t, val, cmap in zip(triangul, values, ['bwr', 'bwr'])] #'bwr' coolwarm


# making the figure
plt.plot([8.5, 8.5], [-0.5, 20.5], color='black', linewidth=2)
plt.plot([12.5, 12.5], [-0.5,20.5], color='black', linewidth=2)

#ax.set_title('')
ax.set_xticks(range(M))
ax.set_xticklabels(['24','25','p26','26','27','28','29','30','31','53','54','55','56','76','77','78','79','80','81','82','83'],
                  ha='center')
ax.set_yticks(range(N))
ax.set_yticklabels(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'], ha='center')
ax.xaxis.tick_top()
ax.yaxis.set_tick_params(pad=10)
ax.set_xlabel('Sites', labelpad=20)
ax.xaxis.set_label_position('top')

ax.spines["bottom"].set_linewidth(2)
ax.spines["top"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)

ax.set_ylabel('Amino Acids', labelpad=10)
ax.invert_yaxis()
ax.margins(x=0, y=0)
ax.set_aspect('equal', 'box')  # square cells

plt.tight_layout()
cbar = plt.colorbar(imgs[0], shrink=0.25, ticks=[0,0.1,0.2])
cbar.ax.set_yticklabels(['0%', '10%', '≥20%']) 


plt.show()
fig.savefig(save)

if __name__ == '__main__':
    main()