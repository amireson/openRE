all: output/mb_dt0-1.csv 

output/mb_dt0-1.csv: Richards.py run.py 
	python -B run.py

.PHONY: clean

clean:
	rm -f output/*.csv
