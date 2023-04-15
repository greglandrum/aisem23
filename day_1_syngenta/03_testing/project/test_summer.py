import pytest
import random
from itertools import permutations
from summer import mysum


# simple testing 1
def test_mysum():
    assert mysum([1, 2, 3]) == 6


# simple testing 2
def test_mysum_tuple():
    assert mysum((1, 2, 3)) == 6


# example of error raising
def test_raises_typeerror():
    with pytest.raises(TypeError):
        mysum([1, 2, '3'])


def test_associative_property():
    # Associative property
    assert mysum([mysum([1, 2]), 3]) == mysum([1, mysum([2, 3])])


# toward property driven testing
def test_commutative_property():
    # Commutative property
    assert mysum([1, 2, 3]) == mysum([3, 1, 2])


def test_additive_identity():
    # Additive Identity property
    assert mysum([1, 2, 3]) == mysum([3, 1, 2, 0])


def test_additive_inverse():
    # Additive Inverse
    assert mysum([1, -1]) == 0


def test_closure_property():
    # Closure property
    x = [1, 2]
    y = mysum(x)
    assert isinstance(y, int)

    x = [1., 2.]
    y = mysum(x)
    assert isinstance(y, float)

    x = [1, 2.]
    y = mysum(x)
    assert isinstance(y, float)

########################################################################################################################
# solutions to generalization




@pytest.fixture(scope="module")
def list_of_int_numbers():
    n = 3  # Adjust this value to control the number of random numbers generated
    numbers = [random.randint(1, 100) for _ in range(n)]
    return numbers


def test_commutative_property_general_1(list_of_int_numbers):
    sums = []

    for perm in permutations(list_of_int_numbers):
        sums.append(mysum(perm))

    assert all(x == sums[0] for x in sums)


def test_associative_property_general_1(list_of_int_numbers):
    associations = []
    for perm in permutations(list_of_int_numbers):
        assoc = mysum([mysum(perm[:2]), perm[2]])
        associations.append(assoc)

    assert all(x == associations[0] for x in associations)


def nested_mysum(numbers):
    if len(numbers) == 1:
        return numbers[0]
    else:
        return mysum([nested_mysum(numbers[:-1]), numbers[-1]])


def test_associative_property_general_2(list_of_int_numbers):
    associations = []

    for perm in permutations(list_of_int_numbers):
        assoc = nested_mysum(perm)
        associations.append(assoc)

    assert all(x == associations[0] for x in associations)
