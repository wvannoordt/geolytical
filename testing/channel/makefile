EXE  := program
MAIN := $(wildcard *.cc)
main:
	${CC_HOST} ${ICONFIG} ${MAIN} -o ${EXE} ${LCONFIG}

clean:
	-rm ${EXE}

test: clean main
	./${EXE}
