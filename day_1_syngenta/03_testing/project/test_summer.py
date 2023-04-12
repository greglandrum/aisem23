import pytest
import random
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


# property drivent testing
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


def test_commutative_property():
    # Commutative property
    assert mysum([1, 2, 3]) == mysum([3, 1, 2])


def test_associative_property():
    # Associative property
    assert mysum([mysum([1, 2]), 3]) == mysum([1, mysum([2, 3])])


def test_additive_identity():
    # Additive Identity property
    assert mysum([1, 2, 3]) == mysum([3, 1, 2, 0])


def test_additive_inverse():
    # Additive Inverse
    assert mysum([1, -1]) == 0
