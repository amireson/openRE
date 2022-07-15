all: output/psi_dt1.csv 

output/psi_dt1.csv: run.py
	python -B run.py

.PHONY: clean

clean:
	rm -f output/*.csv
