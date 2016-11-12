all: 
	$(MAKE) -C src all
	R CMD INSTALL ../lvnm

clean:
	rm src/*.o
	rm src/*.so
