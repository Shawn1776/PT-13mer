13mers: main.cpp
	g++ -g -std=c++11 ranMARS.cpp main.cpp -o main
#run:
#	./main
clean:
	#find . -type f | xargs -n 5 touch # fix "make: Warning: File `makefile' has modification time 48 s in the future"
	#rm -rf $(OBJS)
	rm ./13mers

