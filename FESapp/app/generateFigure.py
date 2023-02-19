import matplotlib.pyplot as plt
import base64
from io import BytesIO


def get_graph():
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png)
    graph = graph.decode('utf-8')
    buffer.close()
    return graph

def get_plot(y, title, xAxisName, yAxisName):
    plt.switch_backend('AGG')
    plt.figure(figsize=(10, 5))
    plt.title(title)
    plt.plot(y)
    plt.xlabel(xAxisName)
    plt.ylabel(yAxisName)
    plt.tight_layout()
    graph = get_graph()
    return graph