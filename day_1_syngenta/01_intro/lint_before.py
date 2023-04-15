import math

# Unused import
import os

# Unused variable
unused_var = 10

def calculate_sum(a, b):
    # Comparison with a literal (should use 'is' instead of '==')
    if a == None or b == None:
        return None

    # Unnecessary use of a comprehension
    numbers = [x for x in range(1, 11)]

    sum = 0
    for number in numbers:
        # Redundant parentheses
        sum += (a * number + b)

    return sum

# No newline at the end of the file
result = calculate_sum(2, 3)
print("Result:", result)