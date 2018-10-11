headerpath=$(PWD)/include
srcpath=$(PWD)/src
libpath=$(PWD)/lib
#Main code
#main.exe: $(srcpath)/main.cpp
#	$(CXX) $(CXXFLAGS) -o ./main.exe $(srcpath)/main.cpp -I$(headerpath)

#Dynamic library linked with th object file
printMessage.so: ./printMessage.o
	$(CXX) $(CXXFLAGS) -shared -o ./printMessage.so ./printMessage.o -I$(headerpath)
	rm ./printMessage.o

printMessage.o: $(srcpath)/printMessage.cpp
	$(CXX) $(CXXFLAGS) -c $(srcpath)/printMessage.cpp -I$(headerpath)
clean:
	rm -f $(libpath)/printMessage.so ./*.o

.PHONY: clean
