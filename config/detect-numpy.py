import os.path
import sys

try:
    import numpy
    numpy_path = numpy.__path__[0]
    numpy_include_path = os.path.join(numpy_path, 'core/include')
    array_header = os.path.join(numpy_include_path, 'numpy/arrayobject.h')
    open(array_header, 'r').read()
    sys.stdout.write('-I%s' % numpy_include_path)
except RuntimError as e:
    sys.stderr.write('Error: %s\nCould not find the numpy C headers. Is the numpy Python module installed properly?\n' % str(e))
    sys.exit(1)
