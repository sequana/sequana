there is a pytest.ini file (empty) that has been added on purpose.
The idea is that at the root level, users can simply type::

    python setup.py test

it reads the setup.cfg and use the default options in particular it sets
--cov=sequana

Now, locally this is an issue. Indeed, once in the ./test directory we generally
want to type::

    pytest test_FILE.py --cov=sequana.FILE

but pytest looks for the setup.cfg and overwrite these options. So, we added
an empty pytest.ini that replaces the values in the setup.cfg in the root
directory.


:reference: https://docs.pytest.org/en/latest/customize.html
