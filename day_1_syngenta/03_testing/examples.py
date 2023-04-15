def mysum(numbers):
    result = sum(numbers)
    return result


x = [1, 2, 3]

# manual testing
y = mysum(x)
print(y)

# automatic testing
assert mysum(x) == 6


def test_mysum():
    assert mysum([1, 2, 3]) == 6  # it should be 6


def test_msum_tuple():
    assert sum((1, 2, 3)) == 6, "Should be 6"


def test_mysum_wrong_results():
    assert mysum([1, 2, 2]) == 6  # it should be 6


def test_mysum_wrong_types():
    assert mysum([1, 2, '3']) == 6  # it should be 6


if __name__ == 'main':
    test_mysum()
    test_msum_tuple()
