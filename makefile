# g++ -o CWSI_test CWSI_test.cpp CWSIFunction.cpp

CC = g++
SOURCES := CWSI_test.cpp CWSIFUnction.cpp
OBJS := $(SOURCES: .cpp=.o)

all: CWSI_test

CWSI_test: $(OBJS)
				$(CC) $(CFLAGS) -o CWSI_test $(OBJS) $(FLAGS) $(LIBS) 

.cpp.o: 
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm -f *o CWSI_test