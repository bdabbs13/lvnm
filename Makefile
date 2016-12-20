all: 
	$(MAKE) -C src all
	R CMD INSTALL ../lvnm

clean:
	rm -f src/*.o
	rm -f src/*.so
	rm -f *~
	rm -f */*~
	rm -f *#
	rm -f */*#
