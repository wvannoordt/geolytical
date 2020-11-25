LIB_NAME := geolytical

ifndef OPTLEVEL
OPTLEVEL := 0
endif

ifndef TESTCLEAN
TESTCLEAN := 0
endif

TESTCLEANTARGET :=
ifeq (${TESTCLEAN}, 1)
TESTCLEANTARGET := clean
endif

BASEIDIR  := $(shell pwd)
SRC_DIR   := ${BASEIDIR}/src
LIB_DIR   := ${BASEIDIR}/lib
OBJ_DIR   := ${BASEIDIR}/obj
HDR_DIR   := ${BASEIDIR}/include

IFLAGS := -I${HDR_DIR}

SRC_FILES := $(wildcard ${SRC_DIR}/*.cpp)

HEADER_FILES := $(wildcard ${SRC_DIR}/*.h)

OBJ_FILES   := $(patsubst ${SRC_DIR}/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))


TARGET := ${LIB_DIR}/lib${LIB_NAME}.a

ifndef CC_HOST
CC_HOST := $(shell which g++)
endif

COMPILE_TIME_OPT := 
COMPILE_TIME_OPT += -DBOUNDS_CHECK=0

HOST_FLAGS := -O${OPTLEVEL} -Wno-unknown-pragmas -g -fPIC -fpermissive -std=c++11

export CC_HOST
export ICONFIG=-I${HDR_DIR}
export LCONFIG=-L${LIB_DIR} -l${LIB_NAME}

.PHONY: final

final: setup ${OBJ_FILES}
	${CC_HOST} -fPIC -shared ${COMPILE_TIME_OPT} ${OBJ_FILES} ${IFLAGS} -o ${TARGET}

${OBJ_FILES}: ${OBJ_DIR}/%.o : ${SRC_DIR}/%.cpp
	${CC_HOST} ${HOST_FLAGS} ${COMPILE_TIME_OPT} ${IFLAGS} -c $< -o $@

setup:
	-rm -r ${HDR_DIR}
	mkdir -p ${LIB_DIR}
	mkdir -p ${OBJ_DIR}
	mkdir -p ${HDR_DIR}
	@for hdr in ${HEADER_FILES} ; do \
		ln -s $${hdr} -t ${HDR_DIR};\
	done

test: ${TESTCLEANTARGET} final
	@for fldr in testing/* ; do \
				echo $${fldr}; \
                ${MAKE} -C $${fldr} -f makefile -s test || exit 1; \
        done
	@echo "${LIB_NAME} passed all tests."

clean:
	-rm -r ${LIB_DIR}
	-rm -r ${OBJ_DIR}
	-rm -r ${HDR_DIR}
	for fldr in testing/* ; do \
	            ${MAKE} -C $${fldr} -f makefile clean ; \
	    done
