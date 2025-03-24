from datetime import datetime
import pandas as pd
import plotly.express as px


df = pd.read_csv('../data/embryo_transmission.csv')


# filter out aneuploidies and failed samples
df = df.loc[(df['embryo_genotype'] != 'failed') & (df['embryo_genotype'] != 'XDXBY') & (df['embryo_genotype'] != 'OY')].copy()


# set marker colours and symbols
colours = ["#636EFA", "#BEBEBE", "#FF8C00", "#000000", "#4682B4"]
#symbols = ['circle', 'x', 'square', 'diamond', 'x']
symbols = ['circle', 'circle', 'circle', 'circle', 'circle']


# change size of marker for adults
marker_size_mapping = {
    'unknown': 5,
    'XBXB': 20,
    'XDXB': 20,
    'XBY': 20,
    'XDY': 20
}

df['marker_size'] = [marker_size_mapping[category] for category in df['genotype']]


# set the angle of the 3d plot
# https://plotly.com/python/3d-camera-controls/?_ga=2.129775182.359524835.1692757213-1789717305.1691966008
x_eye = 2
y_eye = 3
z_eye = 1

fig = px.scatter_3d(df, x='kl2_cq', y='qdal_cq', z='protamine_cq', opacity=1,
                    labels={"kl2_cq": "kl2 Cq", "qdal_cq": "qdal Cq", "protamine_cq": "protamine Cq"},
                    color='genotype', color_discrete_sequence=colours,
                    symbol='genotype', symbol_sequence=symbols, size='marker_size')

fig.update_traces(marker=dict(line=dict(width=1, color='DarkSlateGrey')), selector=dict(mode='markers'))

fig.update_layout(width=1000, height=1000, scene_camera_eye=dict(x=x_eye, y=y_eye, z=z_eye))

# get date for figure name
date = datetime.today().strftime('%Y-%m-%d')

# save figure
# https://plotly.com/python/static-image-export/
# note that the plot itself will still be raster: https://community.plotly.com/t/how-to-export-vector-vectorized-plots/22726
fig.write_image(f"../figures/{date}_embryo_genotyping_3d_scatter.pdf")

print('done')


# fin