import math, sys;

# Naming convention violation: class names should use CamelCase
class my_class:
    # Naming convention violation: use ALL_CAPS for constants
    my_const = 3.14

    # Naming convention violation: method names should use snake_case
    def MyMethod(self,x):return math.sin(x)

# Naming convention violation: function names should use snake_case
def g(x):
    return math.cos(x)

# Naming convention violation: single-letter variable names (avoid using 'I' and 'O')
I = 42

def main():
    x=0; y=sys.argv[0]
    my_obj = my_class()
    for i in range(int(3)):
        x+=my_obj.MyMethod(i)*g(i)
        # Naming convention violation: use snake_case for variable names
        resultVal = x
        print("Result: ", resultVal)

if __name__=='__main__':main()


