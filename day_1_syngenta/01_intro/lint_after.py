import math


# Removed unused import: 'os'

# Removed unused variable: 'unused_var'

def calculate_sum(a, b):
    # Fixed comparison with a literal: use 'is' instead of '=='
    if a is None or b is None:
        return None

    # Simplified the list creation: no need for a comprehension
    numbers = list(range(1, 11))

    sum = 0
    for number in numbers:
        # Removed redundant parentheses
        sum += a * number + b

    return sum


# Added newline at the end of the file
result = calculate_sum(2, 3)
print("Result:", result)
