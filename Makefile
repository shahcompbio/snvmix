SNVMix2: SNVMix2.c SNVMix2.h
	gcc -o SNVMix2 SNVMix2.c -lm

all: SNVMix2

clean:
	rm -f SNVMix2



