include ../bin/rdiff_config.sh

all: example

data: 
	./get_data.sh

example: param


threeexamples: data poisson param nonparam

poisson: data
	./run_example.sh poisson

param: data
	./run_example.sh param

nonparam: data
	./run_example.sh nonparam

clean:
	rm -rf  results-param results-nonparam results-poisson example_data.tar.gz data
