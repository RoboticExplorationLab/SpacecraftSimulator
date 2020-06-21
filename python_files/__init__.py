

from python_test import foo


def test1():
    print('this is test1')
    print('calling foo')
    return foo(14)


def test2():
    print('this is test2')
    return 2


__all__ = [test1, test2, foo]
