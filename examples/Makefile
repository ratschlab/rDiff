include ../bin/rdiff_config.sh

all: example

data: 
	./get_data.sh

example: data
	./run_example.sh param

examples:
	./run_example.sh poisson
	./run_example.sh param
	./run_example.sh nonparam

poisson:
	./run_example.sh poisson

param:
	./run_example.sh param

nonparam:
	./run_example.sh nonparam

clean:
	rm -rf  results-param results-nonparam results-poisson #example_data.tar.gz data
