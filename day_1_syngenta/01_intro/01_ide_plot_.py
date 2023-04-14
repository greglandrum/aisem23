import numpy as np
import matplotlib.pyplot as plt

"""
Command Prompt:
Limitation: Running the code in a command prompt will not display the plot. Instead, it may show a textual representation of the plot, which is not very helpful.
Basic Python Environment (e.g., Jupyter Notebook or IPython):
Advantage: The plot will be displayed inline, allowing for better visualization and exploration.
Drawback: Limited functionality compared to a full-featured IDE.
Integrated Development Environment (IDE) like PyCharm or Visual Studio Code:
Advantage: The IDE provides an interactive environment for running the code and displaying the plot. Users can easily adjust parameters, experiment with different functions, and take advantage of features like code completion, syntax highlighting, and error detection.
Best Practice: Use a powerful IDE for better visualization, debugging, and development experience.
You can use this example in your presentation to demonstrate the benefits of using an IDE for Python development, especially for tasks that involve visualization and interactivity.
"""


def sigmoid(x, a=1, b=0):
    return 1 / (1 + np.exp(-a * (x - b)))


xarr = np.linspace(-10, 10, 100)
yarr = sigmoid(xarr)

plt.plot(xarr, yarr)
plt.xlabel('x')
plt.ylabel('sigmoid(x)')
plt.title('Sigmoid Function')
plt.grid()
plt.show()
