all: result.csv

result.csv: run.py celia_MPM.py
	python -B run.py

.PHONY: clean

clean:
	rm -f result.csv
	rm -f psi.npy
