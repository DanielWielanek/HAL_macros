CC    = g++
LD    = g++
COPTS    = `root-config --cflags` -I/usr/local/root/include -I/home/daniel/Cpp/Makra/HAL_Framework/inst/include -g
LDOPTS    = `root-config --libs` -L//home/daniel/Cpp/Makra/HAL_Framework/inst/lib -lHalAna -lHalCut -lHalData -lHalFeatures -lHalFemto -lHalDbg -lHalOTF -g
#C++ Files
SOURCES =  main.cpp
OBJECTS = $(SOURCES:.cpp=.o)
#Dictionary classes
EXECUTABLE=main.exe
all: $(EXECUTABLE)
$(EXECUTABLE): $(OBJECTS)
	$(LD) -o $@ $^ $(LDOPTS)
#C++ files
.cpp.o:
	$(CC) -o $@ $^ -c $(COPTS)
#Dictionary for ROOT classes
clean:;         @rm -f $(OBJECTS)  $(EXECUTABLE) *.o *.d Dict.cpp Dict.h
