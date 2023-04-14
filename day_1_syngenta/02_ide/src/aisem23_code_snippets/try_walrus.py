from timeit import timeit


def cube(num):
    return num**3


def filter_cubes(num_list, threshold):
    return [cube(x) for x in num_list if cube(x) < threshold]


def filter_cubes_walrus(num_list, threshold):
    return [y for x in num_list if (y:= cube(x)) < threshold]


# try to run this code in the environment with Python 3.7
# and then in the environment with Python 3.9
if __name__ == "__main__":

    num_list = [1, 2, 3, 4, 5, 6, 7, 8]

    result = filter_cubes(num_list, 20)
    print(result)

    result2 = filter_cubes_walrus(num_list, 20)
    print(result2)

    # let's time it and see how long it takes to compute
    time_old = timeit('filter_cubes([1, 2, 3, 4, 5, 6, 7, 8], 2000)',
                      globals=globals())

    time_walrus = timeit('filter_cubes_walrus([1, 2, 3, 4, 5, 6, 7, 8], 2000)',
                         globals=globals())

    print(f"Time for standard way: {time_old:.2f}s")
    print(f"Time for walrus way: {time_walrus:.2f}s")
