import math
import sys


# Naming convention fixed: class names should use CamelCase
class MyClass:
    # Naming convention fixed: use ALL_CAPS for constants
    MY_CONST = 3.14

    # Naming convention fixed: method names should use snake_case
    def my_method(self, x):
        return math.sin(x)


# Naming convention fixed: function names should use snake_case
def cosine_function(x):
    return math.cos(x)


# Naming convention fixed: use more descriptive variable names (avoid using 'I' and 'O')
item_count = 42


def main():
    x = 0
    y = sys.argv[0]
    my_obj = MyClass()

    for i in range(int(3)):
        x += my_obj.my_method(i) * cosine_function(i)
        # Naming convention fixed: use snake_case for variable names
        result_value = x
        print("Result:", result_value)


if __name__ == '__main__':
    main()
