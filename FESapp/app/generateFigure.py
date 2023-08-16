import matplotlib.pyplot as plt
import base64
from io import BytesIO
from cycler import cycler
import matplotlib.lines as mlines
import matplotlib.ticker as plticker




def get_graph():
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png)
    graph = graph.decode('utf-8')
    buffer.close()
    return graph

def get_plot_auto(y, title, xAxisName, yAxisName):
    plt.switch_backend('AGG')
    plt.figure(figsize=(10, 5))
    plt.title(title)
    plt.plot(y)
    plt.xlabel(xAxisName)
    plt.ylabel(yAxisName)
    plt.xlim([0, 120])
    plt.tight_layout()
    graph = get_graph()
    return graph

def get_plot(y, title, xAxisName, yAxisName):
    plt.switch_backend('AGG')
    plt.figure(figsize=(10, 5))
    plt.title(title)
    plt.plot(y)
    plt.xlabel(xAxisName)
    plt.ylabel(yAxisName)
    
    plt.xlim([0, 48])
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(4))

    plt.tight_layout()
    graph = get_graph()
    return graph

#Get plot from multiple data
def get_plot_multiple(y, title, xAxisName, yAxisName):
    plt.switch_backend('AGG')
    plt.figure(figsize=(10, 5))
    plt.title(title)


    blue_star = mlines.Line2D([], [], color='blue', marker='*', linestyle='None',
                          markersize=5, label='Dose 1')
    red_square = mlines.Line2D([], [], color='red', marker='s', linestyle='None',
                            markersize=5, label='Dose 2')
    black_triangle = mlines.Line2D([], [], color='black', marker='^', linestyle='None',
                            markersize=5, label='Dose 3')
    green_circle = mlines.Line2D([], [], color='green', marker='o', linestyle='None',
                            markersize=5, label='Dose 4')
    purple_diamond = mlines.Line2D([], [], color='purple', marker='D', linestyle='None',
                            markersize=5, label='Dose 5')  
    orange_x = mlines.Line2D([], [], color='orange', marker='x', linestyle='None',
                            markersize=5, label='Dose 6')
    #yellow_plus = mlines.Line2D([], [], color='yellow', marker='+', linestyle='None', markersize=10, label='Yellow +')
    #cyan_cross = mlines.Line2D([], [], color='cyan', marker='X', linestyle='None', markersize=10, label='Cyan X')

    
    
    
    

    plt.legend(handles=[blue_star, red_square, purple_diamond, green_circle, black_triangle, orange_x])


    plt.plot(y[0], color='blue', marker='^', linestyle='solid', linewidth=1, markersize=3)
    plt.plot(y[1], color='red', marker='', linestyle='dashdot', linewidth=2, markersize=3)
    plt.plot(y[2], color='black', marker='', linestyle='dashed', linewidth=2, markersize=3)
    plt.plot(y[3], color='green', marker='o', linestyle='solid', linewidth=1, markersize=3)
    plt.plot(y[4], color='purple', marker='D', linestyle='dotted', linewidth=2, markersize=3)
    plt.plot(y[5], color='orange', marker='x', linestyle='solid', linewidth=1, markersize=3)


    plt.xlim([0, 48])
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(4))


    plt.xlabel(xAxisName)
    plt.ylabel(yAxisName)
    plt.tight_layout()
    graph = get_graph()
    return graph