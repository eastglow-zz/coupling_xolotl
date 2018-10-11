headerpath=$(PWD)/include
srcpath=$(PWD)/src
libpath=$(PWD)/lib
#Main code
main.exe: $(srcpath)/main.cpp
	#$(CXX) $(CXXFLAGS) -o ./main.exe $(srcpath)/main.cpp -I$(headerpath)
	g++ $(CXXFLAGS) -o ./main.exe $(srcpath)/main.cpp -I$(headerpath)

#Dynamic library created
triangle.so: ./triangle.o
	#$(CXX) $(CXXFLAGS) -shared -o ./triangle.so ./triangle.o -I$(headerpath)
	$(CXX) $(CXXFLAGS) -shared -o ./triangle.so ./triangle.o -I$(headerpath)
	rm ./triangle.o

triangle.o: $(srcpath)/triangle.cpp
	#$(CXX) $(CXXFLAGS) -c $(srcpath)/triangle.cpp -I$(headerpath)
	g++ $(CXXFLAGS) -c $(srcpath)/triangle.cpp -I$(headerpath)
clean:
	rm -f $(libpath)/triangle.so ./*.o ./main.exe

.PHONY: clean
