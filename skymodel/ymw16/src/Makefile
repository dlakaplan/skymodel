clean:
	rm -f -r build/
	rm -f ymw16_dmdtau.so

clean-c:
	rm -f ymw16_dmdtau.c

clean-all: clean clean-c


.PHONY: build
build: clean
	python setup.py build_ext --inplace

.PHONY: cython-build
cython-build: clean clean-c
	python setup.py build_ext --inplace --use-cython
