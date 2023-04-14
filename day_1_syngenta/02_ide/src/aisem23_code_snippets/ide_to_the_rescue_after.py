import pandas as pd


# in this function we're calculating a square of a given number
def my_square(x):

    square = x * x
    return square


if __name__ == "__main__":

    my_list = [['a', 'b', 'c'], ['AA', 'BB', 'CC']]

    df = pd.DataFrame(my_list, columns=["data_a", "data_b", "data_c"])

    print(df.head(5))

    results = [my_square(x) for x in [1, 2, 3, 4, 5, 6]]

    print(results)

    for item in results:
        w = item + 10
        r = item + 100
        print(r, w)

    for item in results:
        if item == 0:
            print(f"We found zero!!! {item}")
        item += 2
        print(item)
